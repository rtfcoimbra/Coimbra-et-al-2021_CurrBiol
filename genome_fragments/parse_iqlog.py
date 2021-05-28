from sys import argv

infile = argv[1]
is_au = False

with open(infile) as fin:
    for line in fin.readlines():
        if line.startswith('TreeID') and not is_au:
            is_au = True
        elif is_au and line.startswith('Time'):
            is_au = False
            pass
        elif is_au:
            line_split = line.split('\t')
            topology = f"Top{line_split[0]}"
            p_au = line_split[1]
            print(f"{infile}\t{topology}\t{p_au}")
