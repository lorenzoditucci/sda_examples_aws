## headers bug

In file vadd.h, there is a

'''
#define VAL 2
'''

This value is use both in the host and kernel file and it is added to produce the results.

# steps to reproduce the bug
- Build the code for cpu emulation and run it
- change the value of VAL in vadd.h to 0
- Build the code again and run it, the results are now wrong.
