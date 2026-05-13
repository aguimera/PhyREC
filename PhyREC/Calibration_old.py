
import sys
from neo.core import AnalogSignal
import numpy as np
import quantities as pq
from scipy.interpolate import interpolate
from PhyREC.SignalProcess import sliding_window

def CalcVgeff(Sig, Tchar, VgsExp=None, Regim='hole', CalType='interp'):
    """
    Calibrate signal to gate voltage using transfer characteristics.

    Convert measured current signal to effective gate voltage using interpolation
    or linear approximation based on device transfer characteristics.

    Parameters
    ----------
    Sig : neo.core.AnalogSignal
        Input current signal to calibrate.
    Tchar : object
        Device transfer characteristics object with methods GetVgs(), GetIds(),
        GetGM(), GetUd0(), and IsOK attribute.
    VgsExp : float or array-like, optional
        Gate voltage value(s) for linear calibration method.
    Regim : str, optional
        Device regime ('hole' or 'electron', default: 'hole').
    CalType : str, optional
        Calibration method: 'interp' for interpolation or 'linear' for linear
        approximation (default: 'interp').

    Returns
    -------
    neo.core.AnalogSignal
        Calibrated signal in volts with annotations including calibration status,
        bias current, offset current, and device parameters.

    Notes
    -----
    The function adds several annotations to the returned signal:
    - Calibrated: Boolean indicating successful calibration
    - Working: Same as Calibrated
    - IdsOff: Offset current
    - VgsCal: Mean calibrated gate voltage
    - IdsBias: Bias current
    - IsOK: Device transfer characteristics status
    - Iname: Original signal name
    - GM: Conductance in Siemens
    """
    Vgs = Tchar.GetVgs()
    vgs = np.linspace(np.min(Vgs), np.max(Vgs), 10000)

    if Regim == 'hole':
        Inds = np.where(vgs < Tchar.GetUd0())[1]
    else:
        Inds = np.where(vgs > Tchar.GetUd0())[1]

    Ids = Tchar.GetIds(Vgs=vgs[Inds]) * pq.A
    GM = Tchar.GetGM(Vgs=VgsExp) * pq.S

    IdsExp = Tchar.GetIds(Vgs=VgsExp) * pq.A
    IdsOff = np.mean(Sig) - IdsExp
    IdsBias = np.mean(Sig)

    Calibrated = np.array((True,))
    try:
        if CalType == 'interp':
            fgm = interpolate.interp1d(Ids[:, 0], vgs[Inds])
            st = fgm(np.clip(Sig, np.min(Ids), np.max(Ids))) * pq.V
        elif CalType == 'linear':
            st = Sig / GM
        else:
            print('Calibration Not defined')
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape)
        Calibrated = np.array((False,))

    print(str(Sig.name), '-> ',
          'IdsBias', IdsBias,
          'IdsOff', IdsOff,
          'Vgs', np.mean(st),
          Tchar.IsOK)

    annotations = {'Calibrated': Calibrated,
                   'Working': Calibrated,
                   'IdsOff': IdsOff.flatten()[0],
                   'VgsCal': np.mean(st),
                   'IdsBias': IdsBias.flatten()[0],
                   'IsOK': Tchar.IsOK,
                   'Iname': Sig.name,
                   'GM': GM,
                   }

    CalSig = AnalogSignal(st,
                          units='V',
                          t_start=Sig.t_start,
                          sampling_rate=Sig.sampling_rate,
                          name=str(Sig.name),
                          file_origin=Sig.file_origin)

    #    CalSig.annotate(**annotations)
    CalSig.array_annotate(**annotations)

    return CalSig


def CalcVgeff2(Sig, Tchar, VgsExp, Regim='hole', CalType='interp'):
    """
    Calibrate signal to gate voltage using transfer characteristics (alternative).

    Convert measured current signal to effective gate voltage using interpolation
    or linear approximation. This variant flattens arrays and handles units differently
    than CalcVgeff.

    Parameters
    ----------
    Sig : neo.core.AnalogSignal
        Input current signal to calibrate.
    Tchar : object
        Device transfer characteristics object with methods GetVgs(), GetIds(),
        GetGM(), GetUd0(), and IsOK attribute.
    VgsExp : float or array-like
        Gate voltage value(s) for calibration (required, not optional).
    Regim : str, optional
        Device regime ('hole' or 'electron', default: 'hole').
    CalType : str, optional
        Calibration method: 'interp' for interpolation or 'linear' for linear
        approximation (default: 'interp').

    Returns
    -------
    neo.core.AnalogSignal
        Calibrated signal in volts with array annotations including calibration
        status, bias current, offset current, and device parameters.

    Notes
    -----
    Similar to CalcVgeff but uses flattened arrays and stores annotations as
    Python floats rather than quantities objects. Used for compatibility with
    specific analysis workflows.
    """
    Vgs = Tchar.GetVgs()
    vgs = np.linspace(np.min(Vgs), np.max(Vgs), 10000)

    if Regim == 'hole':
        Inds = np.where(vgs < Tchar.GetUd0())[1]
    else:
        Inds = np.where(vgs > Tchar.GetUd0())[1]

    IdsBias = np.array((np.nan,)) * pq.A
    IdsOff = np.array((np.nan,)) * pq.A
    GM = np.nan * pq.S
    try:
        Ids = Tchar.GetIds(Vgs=vgs[Inds]).flatten()
        GM = Tchar.GetGM(Vgs=VgsExp).flatten()

        IdsExp = Tchar.GetIds(Vgs=VgsExp).flatten()
        IdsOff = np.mean(Sig) - IdsExp
        IdsBias = np.mean(Sig)

        Calibrated = np.array((True,))

        if CalType == 'interp':
            fgm = interpolate.interp1d(Ids, vgs[Inds])
            st = fgm(np.clip(Sig, np.min(Ids), np.max(Ids))) * pq.V
        elif CalType == 'linear':
            st = Sig / GM
        else:
            print('Calibration Not defined')
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape) * pq.V
        Calibrated = np.array((False,))

    print(str(Sig.name), '-> ',
          'IdsBias', IdsBias,
          'IdsOff', IdsOff,
          'Vgs', np.mean(st),
          'GM', GM,
          Tchar.IsOK)

    annotations = {'Calibrated': Calibrated[0],
                   'Working': Calibrated[0],
                   'IdsOff': float(IdsOff.flatten()[0].magnitude),
                   'VgsCal': float(np.mean(st).magnitude),
                   'IdsBias': float(IdsBias.flatten()[0].magnitude),
                   'IsOK': Tchar.IsOK,
                   'Iname': Sig.name,
                   'GM': float(GM.magnitude),
                   }

    CalSig = AnalogSignal(st,
                          units='V',
                          t_start=Sig.t_start,
                          sampling_rate=Sig.sampling_rate,
                          name=str(Sig.name),
                          file_origin=Sig.file_origin)

    CalSig.array_annotate(**annotations)
    # CalSig.array_annotate(**annotations)

    return CalSig


def CalcVgeffNoInterp(Sig, Tchar, VgsExp=None, Regim='hole'):
    """
    Calibrate signal to gate voltage using linear approximation only.

    Convert measured current signal to effective gate voltage using conductance-based
    linear calibration. This simplified method does not use interpolation.

    Parameters
    ----------
    Sig : neo.core.AnalogSignal
        Input current signal to calibrate.
    Tchar : object
        Device transfer characteristics object with method GetGM() and IsOK attribute.
    VgsExp : float or array-like, optional
        Gate voltage value(s) for conductance calculation.
    Regim : str, optional
        Device regime ('hole' or 'electron', default: 'hole'). Not currently used.

    Returns
    -------
    neo.core.AnalogSignal
        Calibrated signal in volts with array annotations including calibration
        status, mean gate voltage, and device status.

    Notes
    -----
    This is a simplified calibration that only uses conductance (GM) for linear
    scaling. No interpolation is performed on transfer characteristics.
    """
    gm = Tchar.GetGM(Vgs=VgsExp)

    Calibrated = np.array((True,))
    try:
        st = Sig.magnitude / gm
    except:
        print(Sig.name, "Calibration error:", sys.exc_info()[0])
        st = np.zeros(Sig.shape)
        Calibrated = np.array((False,))

    print(str(Sig.name), '-> ', 'GM', gm, 'Vgs', np.mean(st), Tchar.IsOK)
    annotations = {'Calibrated': Calibrated,
                   'Working': Calibrated,
                   # 'IdsOff': IdsOff.flatten(),
                   'VgsCal': np.array((np.mean(st),)),
                   'IsOK': np.array((Tchar.IsOK,)),
                   'Iname': np.array((Sig.name,)),
                   }

    CalSig = AnalogSignal(st * pq.V,
                          units='V',
                          t_start=Sig.t_start,
                          sampling_rate=Sig.sampling_rate,
                          name=str(Sig.name),
                          file_origin=Sig.file_origin)

    #    CalSig.annotate(**annotations)
    CalSig.array_annotate(**annotations)

    return CalSig


def NoiseBlanking(sig, NoiseLimitRMS=20 * pq.uV, MinNoiseFreeTime=5 * pq.s, timewidth=1 * pq.s, Value=0):
    """
    Remove high-noise sections from signal by blanking with a specified value.

    Identify periods where signal RMS exceeds a threshold and replace them with
    a constant value. Maintains breaks between noise blocks if they exceed minimum
    duration.

    Parameters
    ----------
    sig : neo.core.AnalogSignal
        Input signal to process.
    NoiseLimitRMS : quantities.Quantity, optional
        RMS threshold for noise detection (default: 20 µV).
    MinNoiseFreeTime : quantities.Quantity, optional
        Minimum time between noise blocks required to split into separate blocks
        (default: 5 s).
    timewidth : quantities.Quantity, optional
        Time window for sliding RMS calculation (default: 1 s).
    Value : float, optional
        Value to replace noise sections with (default: 0).

    Returns
    -------
    neo.core.AnalogSignal
        Signal with noise regions blanked with the specified value.

    Notes
    -----
    The function uses sliding window RMS analysis to identify noisy periods.
    Consecutive noise blocks separated by less than MinNoiseFreeTime are merged.
    """
    sRMS = sliding_window(sig, timewidth)

    noisets = sRMS.times[np.where(sRMS > NoiseLimitRMS)[0]]

    if len(noisets) == 0:
        return sig

    inds = np.where(np.diff(noisets) > MinNoiseFreeTime)[0]

    NoiseBlocks = []
    NoiseBlocks.append((noisets[0], noisets[inds[0]]))
    for ic, ind in enumerate(inds[:-1]):
        NoiseBlocks.append((noisets[ind + 1], noisets[inds[ic + 1]]))
    NoiseBlocks.append((noisets[inds[-1] + 1], noisets[-1]))

    for ti, te in NoiseBlocks:
        i1 = sig.time_index(ti)
        i2 = sig.time_index(te)
        sig[i1:i2] = np.ones((i2 - i1, 1)) * Value * sig.units

    return sig
