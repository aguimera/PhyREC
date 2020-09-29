# -*- coding: utf-8 -*-
# Copyright © 2016, German Neuroinformatics Node (G-Node)
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted under the terms of the BSD License. See
# LICENSE file in the root of the Project.

from enum import Enum


class LinkType(Enum):
    # Values of enum are also part of NIX format
    Tagged = "tagged"
    Untagged = "untagged"
    Indexed = "indexed"
