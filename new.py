import argparse
parser = argparse.ArgumentParser(description="function: return a * b=? ")
parser.add_argument('-F', "--first", required=True, type=int, help="input a number")
parser.add_argument('-S', "--second", required=True, type=int, help="input a number")
args = parser.parse_args()

F=args.first
S=args.second
def my_abs(F,S):
    print  F*S
my_abs(F,S)