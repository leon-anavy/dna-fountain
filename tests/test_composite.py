#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
# Licensed under the MIT License.
# Filename: test_composite
# Project: dna-fountain-clean
# Author: Leon Anavy
# Email: anavy (at) technion (dot) ac (dot) il
# Created: 19/02/2019 14:32

from utils import *
import random, sys

alphabet_file = '/data/anavy/Storage/composite_fountain/simulate_phi5/alphabet_phi5.txt'

alphabet = read_composite_alphabet(alphabet_file)
bin_block_size,output_block_size,comp_encoder,comp_decoder = create_composite_encoder(alphabet,23,4)


# random intarray
for _ in range(100):
	data_int = [random.randint(0,255) for _ in range(92)]
	data_bin = int_array_to_bin(data_int)
	data_DNA = [comp_encoder(data_bin[i:i + bin_block_size]) for i in range(0, len(data_bin), bin_block_size)]
	data_bin1 = ''.join([comp_decoder(d) for d in data_DNA])
	print data_bin1 == data_bin

# random binary string
for _ in range(100):
	data_bin = ''.join([str(random.randint(0,1)) for _ in range(23)])
	data_DNA = [comp_encoder(data_bin[i:i + bin_block_size]) for i in range(0, len(data_bin), bin_block_size)]
	data_bin1 = [comp_decoder(d) for d in data_DNA][0]
	print data_bin1 == data_bin
