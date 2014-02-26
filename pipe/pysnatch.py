# -*- coding: utf-8 -*-
"""
pysnatch.py is licensed under The MIT License (MIT)

Copyright (c) 2014 Eric J. Mott

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Created on Fri Feb  7 16:49:17 2014
@author: Mott
"""
from numpy import matrix
from xlrd import open_workbook

def xlsnatch(file_name_as_string,sheet_index,starting_row,starting_col,num_of_rows,num_of_cols): 
    """
    Imports data from a specified sheet from an MS Excel Workbook as a num_of_rows x num_of_cols matrix.
    
    Dependencies: NumPy and xlrd
    """
    
    sheet = open_workbook(file_name_as_string).sheet_by_index(sheet_index)
    
    if num_of_rows > sheet.nrows - starting_row + 1 | num_of_cols > sheet.ncols - starting_col + 1:
        print 'You requested more rows or columns than there are in the sheet!'
    
    return matrix([[(sheet.cell(row_index,col_index).value if sheet.cell(row_index,col_index).value !='' else 0) for col_index in range(starting_col-1,starting_col-1+num_of_cols)] for row_index in range(starting_row-1,starting_row-1+num_of_rows)])
    
