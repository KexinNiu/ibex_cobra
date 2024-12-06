###########################################
#########Evaluation of the model##########
###########################################

# with input of model name and the seed, the function will return the evaluation result of the model
# combined with the recon gapfill or not


import argparse
parser = argparse.ArgumentParser(description='Precentage to remove reactions')
parser.add_argument('--compare', default = True, help='if compare with the RECON gapfill')
parser.add_argument('--model', default = 'iJO1366', help='model name')
parser.add_argument('--prec', default = 0.05, help='percentage of reactions to remove')
parser.add_argument('--seednum', default = 1234, help='seed for random')

args = parser.parse_args()
print('args:',args.prec,args.seednum)