#!/usr/bin/env python

"""
Auxiliary functions to find next simulation id and next experiment in a given directory
"""

import os

def find_next_available_file(word, path):
    next_available = 0
    for entry in os.scandir(path):
        if entry.is_file() and word in entry.name:
            current_str = entry.name[-6:]
            start = 0
        
            length = len(current_str)
            for i in range(length):
                if current_str[i] == '0':
                    start = start + 1
                else:
                    break
           
            current_int = int(current_str[start:]) if len(current_str[start:]) > 0 else 0
            if current_int >= next_available:
                next_available = current_int + 1
    os.scandir(path).close()
    return word + str(next_available).zfill(6)

def find_next_available_dir(word, path):
    next_available = 0
    for entry in os.scandir(path):
        if entry.is_dir() and word in entry.name:
            current_str = entry.name[-6:]
            start = 0
        
            length = len(current_str)
            for i in range(length):
                if current_str[i] == '0':
                    start = start + 1
                else:
                    break
           
            current_int = int(current_str[start:]) if len(current_str[start:]) > 0 else 0
            if current_int >= next_available:
                next_available = current_int + 1
    os.scandir(path).close()
    return word + str(next_available).zfill(6)