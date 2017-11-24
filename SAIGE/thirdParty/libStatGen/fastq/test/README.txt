Lines 1 - 24 - test that all valid quality string characters are accepted & tests multiple line Raw Sequence and Quality Strings.

Sequence Identifier Line Validates:
    * Line 25: line is at least 2 characters long ('@' and at least 1 for the sequence identifier)
    * Line 29: line starts with an '@'
    * Line 33 & 37: no space between the '@' & the sequence identifier (which must be at least 1 character)
    * Line 41: sequence identifier is unique within the file

Raw Sequence Line Validates:
    * Line 46 & 47: every character is in ACTGNactgn0123.
    * Line 51: the raw sequence after it is completely read is at least a configurable minimum length
    * Line 56 & 57: assumes all lines are part of the raw sequence until a line begins with a '+' or the end of the file is reached

Plus Line Validates:
    * Line 88: sequence identifier on + line does not match the one on the @ line.
    * Line 91: that this line exists for each sequence

Quality Line Validates:
    * Line 63 & 64: each character is > ascii 32
    * Line 70: assumes all lines are part of the quality string until the total length of quality characters is >= the raw sequence length or the end of the file is reached
    * Line 77: length of the quality string equals the length of the raw sequence


