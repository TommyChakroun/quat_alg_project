import ast

def read_line(string):
    """
    Parse a line:
    - If it's an integer, return int.
    - If it's a list of lists, parse manually and return as list of list of QQ.
    """
    s = string.strip()

    # Case: integer
    if s[0] != "[":
        return QQ(s)

    import re
    try:
        # Match all lists inside the outer list using regex
        row_strings = re.findall(r'\[([^\[\]]+)\]', s)
        matrix = []
        for row_str in row_strings:
            entries = [QQ(entry.strip()) for entry in row_str.split(',')]
            matrix.append(entries)
        return matrix
    except Exception as e:
        raise ValueError(f"Could not parse matrix string: {e}")


def read_one_block(lines):
    """
    Reads a block starting at start_index.
    The block alternates: int, list_of_lists, ..., and may end in either one.
    Returns (block_data, next_index).
    """
    return [read_line(line) for line in lines if line.strip()]


def read_file(size_block,filename):
    """
    Return list_of_blocks the list of each block
    """
    list_blocks = []
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]  # skip empty lines
    
    i = 0
    while i + size_block+1 < len(lines):
        block = []
        for k in range(size_block):
            block.append(lines[i+k])
        list_blocks.append(read_one_block(block))
        i = i+size_block
    return list_blocks