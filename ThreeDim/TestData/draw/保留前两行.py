import os

fileInPath = r"step_two_fliter_10w02061935.dat"
fileOutPath = r"New_step_two_10w02060859.dat"


def main():
	fileIn = open(fileInPath, 'r')
	fileLines = fileIn.readlines()
	newLines = [x.strip('\n').split('\t')[6,12] for x in fileLines]
	outLines = ['\t'.join(x)+'\n' for x in newLines]
	fileOut = open(fileOutPath, 'w')
	fileOut.writelines(outLines)

if __name__ == '__main__':
	main()
