#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import neo
import numpy as np
import quantities as pq
import os

import xml.etree.ElementTree as ET
        
def GetNeoItem (List, ItemName):
    ItemNames = {}
    for i, item in enumerate(List):
        name = item.name
        if item.name is None:
            print('Item ', i, ' No name' )
            continue           
        if name in ItemNames.keys():
            print('Item ', i, name, 'duplicated name' )
        ItemNames.update({name: item})
    
    if ItemName in ItemNames.keys():
        return ItemNames[ItemName]
    else:
        print (ItemName, 'No Found')
        return None

        
def dict_to_xml(tag, d):
    '''
    Turn a simple dict of key/value pairs into XML
    '''
    elem = ET.Element(tag)
    for key, val in d.items():
        child = ET.Element(key)
        child.text = str(val)
        child.tail = '\n'
        elem.append(child)
    elem.tail = '\n'
    return elem


class NeoTrain(neo.SpikeTrain):

    def CheckTime(self, Time):
        if Time is None:
            return (self.t_start, self.t_stop)

        if len(Time) == 1:
            Time = (Time[0], Time[0] + self.sampling_period)

        if Time[0] is None or Time[0] < self.t_start:
            Tstart = self.t_start
        else:
            Tstart = Time[0]

        if Time[1] is None or Time[1] > self.t_stop:
            Tstop = self.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def GetSignal(self, Time, Units='s'):
        time = self.CheckTime(Time)
        sl = self.time_slice(time[0], time[1])
        return sl.rescale(Units)


class NeoSignal(neo.AnalogSignal):
    ProcessChain = None
    ProcessChainTime = None

    def CheckTime(self, Time):
        if Time is None:
            return (self.t_start, self.t_stop)

        if len(Time) == 1:
            Time = (Time[0], Time[0] + self.sampling_period)

        if Time[0] is None or Time[0] < self.t_start:
            Tstart = self.t_start
        else:
            Tstart = Time[0]

        if Time[1] is None or Time[1] > self.t_stop:
            Tstop = self.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def GetSignal(self, Time, Units=None):        

        if self.ProcessChain is None:
            time = self.CheckTime(Time)
            sl = self.time_slice(time[0], time[1])
            if Units is not None:
                sl = sl.rescale(Units)
            return sl
        else:
            ProcTime = self.CheckTime(self.ProcessChainTime)
            sl = self.time_slice(ProcTime[0], ProcTime[1])
            for Proc in self.ProcessChain:
                sl = Proc['function'](sl, **Proc['args'])

            sl.__class__ = NeoSignal
            time = sl.CheckTime(Time)
            sl = sl.time_slice(time[0], time[1])
            if Units is not None:
                sl = sl.rescale(Units)
            return sl

#    def AppendSignal(self, Vect):
#        v_old = np.array(self.Signal)
#        v_new = np.vstack((v_old, Vect))*self.Signal.units
#
#        self.Signal = self.Signal.duplicate_with_new_array(signal=v_new)


class NeoSegment():
    def __init__(self, RecordFile=None, Seg=None):
        self.SigNames = {}
        self.EventNames = {}
        if RecordFile is None:
            if Seg is None:
                self.Seg = neo.Segment()
            else:
                self.Seg = Seg
                self.UpdateEventDict()
                self.UpdateSignalsDict()
            return

        ftype = RecordFile.split('.')[-1]
        if ftype == 'h5':
            self.RecFile = neo.io.NixIO(filename=RecordFile, mode='ro')
            Block = self.RecFile.read_block()
        elif ftype == 'smr':
            self.RecFile = neo.io.Spike2IO(filename=RecordFile)
            Block = self.RecFile.read()[0]

        self.Seg = Block.segments[0]

        self.UpdateSignalsDict()
        self.UpdateEventDict()

    def ExportNeuroscope(self, FileName, Range, Bits, Units,
                         ChNames, NeuroScopeMap,
                         ProcessChain=None):     
       
        LSB = (Range / 2)/(2**(Bits-1))
        
        ExpDat = np.ndarray((self.GetSignal(ChNames[0]).size,
                             len(ChNames)),
                            dtype=np.int16)
        
        for ich, chn in ChNames.items():
            sig = self.GetSignal(chn)
            if ProcessChain is not None:
                sig.ProcessChain = ProcessChain
            dat = sig.GetSignal(None, Units=Units)
            dexp = np.array(dat)/LSB
            dexp = np.array(dexp).astype(np.int16)            
            ExpDat[:, ich] = dexp.flatten()
        
        ExpDat.astype(np.int16).tofile(FileName + '.dat')

        Etop = ET.Element('parameters')
        Etop.set('version', '1.0')
        Etop.set('creator', 'neuroscope-2.0.0')
        Acq = {'nBits' : Bits,
               'nChannels': ExpDat.shape[1],
               'samplingRate': str(sig.sampling_rate.magnitude),
               'voltageRange': Range,
               'amplification': str(1),
               'offset': str(0)}
        
        Eacq = dict_to_xml('acquisitionSystem', Acq)
        
        Etop.append(Eacq)

        Afp = {'lfpSamplingRate': str(1250),}
        Eafp = dict_to_xml('fieldPotentials',Afp)
        Etop.append(Eafp)
        
        
        Echg = ET.Element('channelGroups')
        Echg.tail = '\n'
        for g in NeuroScopeMap:
            Eg = ET.Element('group')
            Eg.tail = '\n'
            for chn, ich in zip(g[0], g[1]):
                Ech = ET.Element('channel',
                                 attrib={'skip': str(0),
                                         'name': chn})
                Ech.text = str(ich)
                Ech.tail = '\n'
                Eg.append(Ech)
            Echg.append(Eg)
        
        Eanadesc = ET.Element('anatomicalDescription')
        Eanadesc.tail = '\n'
        Eanadesc.append(Echg)
        Etop.append(Eanadesc)
        
        Espike = ET.Element('spikeDetection')
        Espike.tail = '\n'
        Etop.append(Espike)
        
        
        TopTree = ET.ElementTree(element=Etop)

#        print TopTree
#        print ET.dump(TopTree)
        TopTree.write(FileName + '.xml',
                       encoding="utf-8", xml_declaration=True)




    def SaveRecord(self, FileName, OverWrite=True):
        if os.path.isfile(FileName):
            if OverWrite:
                os.remove(FileName)
            else:
                print ('Warning File Exsist')

        if FileName.endswith('.h5'):
            out_f = neo.io.NixIO(filename=FileName)
        elif FileName.endswith('.mat'):
            out_f = neo.io.NeoMatlabIO(filename=FileName)
        else:
            return

        out_bl = neo.Block(name='NewBlock')
        out_bl.segments.append(self.Seg)
        out_f.write_block(out_bl)
        if FileName.endswith('.h5'):
            out_f.close()

    def UpdateSignalsDict(self):
        self.SigNames = {}
        for i, sig in enumerate(self.Seg.analogsignals):
            if sig.name is None:
                name = str(i)
                if sig.annotations is not None:
                    if 'nix_name' in sig.annotations.keys():
                        name = sig.annotations['nix_name']
            else:
                name = str(sig.name)

            self.SigNames.update({name: i})
        self.signames = self.SigNames.keys()

    def UpdateEventDict(self):
        self.EventNames = {}
        for i, eve in enumerate(self.Seg.events):
            if eve.name is None:
                try:
                    name = eve.annotations['title']
                except:
                    print ('Event found no name ', i)
                    name = str(i)
                self.EventNames.update({name: i})
            else:
                self.EventNames.update({eve.name: i})

    def GetEventTimes(self, EventName, Time=None):
        eve = self.Seg.events[self.EventNames[EventName]].times
        if Time:
            events = eve[np.where((eve > Time[0]) & (eve < Time[1]))]
        else:
            events = eve
        return events

    def GetTstart(self, ChName):
        return self.Seg.analogsignals[self.SigNames[ChName]].t_start

    def SetTstart(self, ChName, Tstart):
        self.Seg.analogsignals[self.SigNames[ChName]].t_start = Tstart

    def SetSignal(self, ChName, Sig):
        self.Seg.analogsignals[self.SigNames[ChName]] = Sig

    def GetSignal(self, ChName):
        sig = self.Seg.analogsignals[self.SigNames[ChName]]
        sig.__class__ = NeoSignal
        return sig

    def Signals(self):
        for s in self.Seg.analogsignals:
            s.__class__ = NeoSignal
        return self.Seg.analogsignals

    def AddEvent(self, Times, Name):
        eve = neo.Event(times=Times,
                        units=pq.s,
                        name=Name)

        self.Seg.events.append(eve)
        self.UpdateEventDict()

    def AddSignal(self, Sig):
        self.Seg.analogsignals.append(Sig)
        self.UpdateSignalsDict()

    def AppendSignal(self, ChName, Vect):
        sig = self.Seg.analogsignals[self.SigNames[ChName]]

        v_old = np.array(sig)
        v_new = np.vstack((v_old, Vect))

        sig_new = sig.duplicate_with_new_array(signal=v_new*sig.units)

        self.SetSignal(ChName, sig_new)


#def ReadMCSFile(McsFile, OutSeg=None, SigNamePrefix=''):
#    import McsPy.McsData as McsData
#
#    Dat = McsData.RawData(McsFile)
#    Rec = Dat.recordings[0]
#    NSamps = Rec.duration
#
#    if OutSeg is None:
#        OutSeg = NeoSegment()
#
#    for AnaStrn, AnaStr in Rec.analog_streams.iteritems():
#        if len(AnaStr.channel_infos) == 1:
#            continue
#
#        for Chn, Chinfo in AnaStr.channel_infos.iteritems():
#            print ('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
#            ChName = str(SigNamePrefix + Chinfo.label)
#            print (ChName)
#
#            Fs = Chinfo.sampling_frequency
#            Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)
#            sig = neo.AnalogSignal(pq.Quantity(Var, Chinfo.info['Unit']),
#                                   t_start=0*pq.s,
#                                   sampling_rate=Fs.magnitude*pq.Hz,
#                                   name=ChName)
#
#            OutSeg.AddSignal(sig)
#    return OutSeg
