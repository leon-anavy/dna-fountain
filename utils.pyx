"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from string import maketrans   # Required to call maketrans function.
import struct
import random
import os
import numpy as np
import argparse
import math, itertools

intab = "0123"
outtab = "ACGT"

trantab = maketrans(intab, outtab)
revtab = maketrans(outtab, intab)


def charN(str, N):
    if N < len(str):
        return str[N]
    return 'X'
    
def xor(str1, str2, length):
    return ''.join(chr(ord(str1[i]) ^ ord(str2[i])) for i in xrange(length))


def xor_np(ord_array1, ord_array2):
    return np.bitwise_xor(ord_array1, ord_array2)

def xor_ord(ord_array1, ord_array2, length):
    #ord_array_res = [None] * length
    for i in xrange(length):
        #ord_array_res[i] = ord_array1[i] ^ ord_array2[i]
        ord_array1[i] = ord_array1[i] ^ ord_array2[i]
    return ord_array1#_res


def xor_bin(str1, str2):
    length = max(len(str1),len(str2))
    bytes = ''
    for i in xrange(length):
         byte_xor_dec = ord(charN(str1,i)) ^ ord(charN(str2,i))
         byte_xor_bin = "{0:08b}".format(byte_xor_dec)
         bytes = ''.join([bytes, byte_xor_bin])
    return bytes



def bin_to_dna(bin_str):
    s = ''.join(str(int(bin_str[t:t+2],2)) for t in xrange(0, len(bin_str),2)) #convert binary 2-tuple to 0,1,2,3
    return s.translate(trantab)
 

def byte_to_dna(s):
    #convert byte data (\x01 \x02) to DNA data: ACTC
    bin_data = ''.join('{0:08b}'.format(ord(s[t])) for t in xrange(0,len(s)))
    return bin_to_dna(bin_data)


def byte_from_bin(s):
    #convert string like 01010101 to string of bytes like \x01 \x02 \x03
    return ''.join(chr(int(s[t:t+8],2)) for t in xrange(0, len(s), 8))

def byte_to_int_array(s):
    #convert a strong like \x01\x02 to [1,2,]
    a = list()
    for t in xrange(0, len(s)):
        a.append(ord(s[t]))
    return a

def int_array_to_bin(s):
    return ''.join('{0:08b}'.format(element) for element in s) #convert to a long sring of binary values

def int_to_dna(a):
    #a is an array of integers between 0-255.
    #returns ACGGTC
    bin_data = ''.join('{0:08b}'.format(element) for element in a) #convert to a long sring of binary values
    s = ''.join(str(int(bin_data[t:t+2],2)) for t in xrange(0, len(bin_data),2)) #convert binary array to a string of 0,1,2,3
    return s.translate(trantab)


def int_to_four(a):
    #a is an array of integers between 0-255.
    #returns 0112322102
    bin_data = ''.join('{0:08b}'.format(element) for element in a) #convert to a long sring of binary values
    return ''.join(str(int(bin_data[t:t+2],2)) for t in xrange(0, len(bin_data),2)) #convert binary array to a string of 0,1,2,3

def four_to_dna(s):
    return s.translate(trantab)

def dna_to_byte(dna_str):
    #convert a string like ACTCA to a string of bytes like \x01 \x02
    num = dna_str.translate(revtab)
    s = ''.join('{0:02b}'.format(int(num[t])) for t in xrange(0, len(num),1))
    data = ''.join(chr(int(s[t:t+8],2)) for t in xrange(0, len(s), 8))

    return data

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    num = dna_str.translate(revtab)
    s = ''.join('{0:02b}'.format(int(num[t])) for t in xrange(0, len(num),1))
    data = [int(s[t:t+8],2) for t in xrange(0,len(s), 8)]

    return data

def comp_to_int_array(dna_str,comp):
    if dna_str is None:
        return  None
    #convert a composite dna word to an array of ints like [10, 2, 4]
    alphabet = comp['alphabet']
    rev_alphabet = {v:k for k,v in alphabet.iteritems()}
    comp_encoder = comp['encoder']
    BC_bases = comp['BC_bases']
    std_DNA = dna_str[:BC_bases]
    comp_DNA = dna_str[BC_bases:]
    std_int_ar = dna_to_int_array(std_DNA)
    comp_block_size = len(comp_encoder.values()[0])
    comp_encoder_rev = {v:k for k,v in comp_encoder.iteritems()}
    comp_DNA_ints = [rev_alphabet[l] for l in comp_DNA]
    comp_bin = ''.join(comp_encoder_rev[tuple(comp_DNA_ints[i:i + comp_block_size])] for i in range(0, len(comp_DNA_ints), comp_block_size))
    bytes = [byte_from_bin(comp_bin[i:i + 8]) for i in range(0, len(comp_bin), 8)]
    comp_int_ar = byte_to_int_array(bytes)
    data = std_int_ar + comp_int_ar
    return data


def split_header(data_str, header_bytes):
    data = data_str[header_bytes:]
    header_raw = data_str[:header_bytes]
    header_binary = ''.join('{0:08b}'.format(ord(header_raw[t])) for t in xrange(0,header_bytes))
    
    

    header = 0
    for t in xrange(0, len(header_binary)):

        header = header << 1
        header +=  int(header_binary[t])
    
    return (header, data)

def prepare(max_repeat):
    global As, Cs, Gs, Ts
    As = '0' * (max_repeat+1)
    Cs = '1' * (max_repeat+1)
    Gs = '2' * (max_repeat+1)
    Ts = '3' * (max_repeat+1)

def screen_repeat_composite(drop,comp):
    # This is "fake" translation to DNA only to get the BC
    dna = drop.toDNA(comp=None)
    bc = dna[0:comp['BC_bases']]
    if not screen_repeat_dna(bc,4,0.2):
        return 0
    return 1

def screen_repeat(drop, max_repeat, gc_dev):

    dna = drop.toDNA()
    return screen_repeat_dna(dna, max_repeat, gc_dev)
    
def screen_repeat_dna(dna, max_repeat, gc_dev):
    
    if As in dna or Cs in dna or Gs in dna or Ts in dna: 
        return 0

    gc = dna.count("1") + dna.count("2")  
    gc = gc/(len(dna)+0.0)

    if (gc < 0.5 - gc_dev) or (gc > 0.5 + gc_dev):
        return 0
    return 1


def hamming_distance(x,y):

    d = 0
    for t in xrange(0, max(len(x), len(y))):
        if(x[t] != y[t]):
            d += 1
    return d

def restricted_float(x):
    #helper function from http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def read_composite_alphabet(alphabet_file):
    alphabet = {i:l for i,l in zip(range(4),'ACGT')}
    with open(alphabet_file) as inf:
        for idx, l in enumerate(inf):
            alphabet[idx+4] = l.rstrip()
    return alphabet

def create_composite_encoder(alphabet,max_block,oligo_len):
    alphabet_size = len(alphabet)
    get_bin = lambda x, n: format(x, 'b').zfill(n)
    # calcualte best bloxk sizes
    best_binary_b = 0
    best_output_b = 0
    best_encoded_bits = 0
    for binary_b in range(1, max_block):
        binary_words = 2 ** binary_b
        output_b = int(math.ceil(math.log(binary_words, alphabet_size)))
        output_blocks = int(math.floor(oligo_len / output_b))
        encoded_bits = output_blocks * binary_b
        if encoded_bits > best_encoded_bits:
            best_binary_b = binary_b
            best_output_b = output_b
            best_encoded_bits = encoded_bits

    binary_block_size = best_binary_b
    output_block_size = best_output_b
    if output_block_size > 1:
        output_blocks = list(itertools.product(alphabet.keys(), repeat=output_block_size))
    else:
        output_blocks = alphabet.keys()
    encoding_dict = {get_bin(v, binary_block_size): output_blocks[v] for v in
                          range(2 ** binary_block_size)}
    return encoding_dict





