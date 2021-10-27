import pandas as pd
import argparse
import time
import os
import numpy as np
from scipy.special import softmax
import pickle

np.set_printoptions(precision=4, linewidth=200)

def title():
    print('**********************************************************************')
    print('* SparsePro for efficient genome-wide fine-mapping                   *')
    print('* Version 1.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca                      *')
    print('**********************************************************************')
    print()

def get_XX_XtX_ytX(LD,beta,se,var_Y):
    '''get sufficient statistics from summary statistics'''
    XX = var_Y/(se**2)
    XtX = LD * var_Y / (np.dot(se.reshape(-1,1),se.reshape(1,-1)))
    ytX = XX * beta
    return XX, XtX, ytX

#unstandardized HESS extended from Shi et al.,2016
def get_HESS_h2_SS(XtX,XX,LD,beta,se,N,var_Y,LDthres=0.1):
    '''calculate local heritabilities'''
    idx_retain = []
    idx_exclude = [i for i in range(len(beta))]
    zscore=np.abs(beta/se)
    while len(idx_exclude)>0:
        maxid = idx_exclude[np.argmax(zscore[idx_exclude])]
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(LD[maxid,:]>LDthres)[0]]
    Indidx = np.sort(idx_retain)
    #obtain independent signals
    P = len(Indidx)
    XtX_id = XtX[np.ix_(Indidx,Indidx)]
    R_inv = np.linalg.inv(XtX_id)
    vec_id = XX[Indidx] * beta[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(),R_inv),vec_id)-var_Y*P)/(var_Y*(N-P))
    
    if h2_hess<0.0001:
        h2_hess = 0.0001
    if h2_hess>0.9:
        h2_hess = 0.9
    return h2_hess

class SparsePro(object):
    
    def __init__(self,P,K,XX,var_Y,h2):
        '''initialize and set hyperparameters'''
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p,self.k))
        self.beta_mu = np.zeros((self.p,self.k))
        self.beta_prior_tau = np.tile((1.0 / var_Y * h2 * np.array([k+1 for k in range(self.k)])),(self.p,1))
        self.y_tau = 1.0 / (var_Y * (1-h2))
        self.prior_pi = np.ones((self.p,)) * (1/self.p)
        self.beta_post_tau = np.tile(XX.reshape(-1,1),(1,self.k)) * self.y_tau + self.beta_prior_tau
        
    def infer_q_beta(self,XX,ytX,XtX,LD):
        '''perform variational updates'''
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:,idxall] * self.beta_mu[:,idxall]).sum(axis=1)
            self.beta_mu[:,k] = (ytX-np.dot(beta_all_k, XtX))/self.beta_post_tau[:,k] * self.y_tau
            u = -0.5*np.log(self.beta_post_tau[:,k]) + np.log(self.prior_pi.transpose()) + 0.5 * self.beta_mu[:,k]**2 * self.beta_post_tau[:,k]
            self.gamma[:,k] = softmax(u)
            maxid = np.argmax(u)
            self.gamma[abs(LD[maxid])<0.05,k]= 0.0

    def get_elbo(self):
        
        beta_all = (self.gamma * self.beta_mu).sum(axis=1)
        ll1 = self.y_tau * np.dot(beta_all,ytX)
        ll2 = - 0.5 * self.y_tau * ((((self.gamma * self.beta_mu**2).sum(axis=1) * XX).sum()))
        W = self.gamma * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(),XtX),W)
        ll3 = - 0.5 * self.y_tau * ( WtRW.sum() - np.diag(WtRW).sum())
        ll = ll1 + ll2 + ll3
        betaterm1 = -0.5 * (self.beta_prior_tau * self.gamma * (self.beta_mu**2)).sum()
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1,1),(1,self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma!=0] * np.log(self.gamma[self.gamma!=0])).sum()
        mkl = betaterm1 + gammaterm1 - gammaterm2
        elbo = ll + mkl
        
        return ll,mkl,elbo
       
    def get_PIP(self):
        
        return np.max((self.gamma),axis=1)
        
    def update_pi(self, new_pi):
        
        self.prior_pi = new_pi

    def get_effect_dict(self):
        
        numidx = (self.gamma>0.1).sum(axis=0)
        matidx = np.argsort(-self.gamma, axis=0)
        return {i:matidx[0:numidx[i],i].tolist() for i in range(self.k) if numidx[i]>0}
    
    def get_effect_num_dict(self):
        
        gamma = np.round(self.gamma,4)
        beta_mu = np.round(self.beta_mu,4)
        effect = self.get_effect_dict()
        eff_gamma = {i:gamma[effect[i],i].tolist() for i in effect}
        eff_mu = {i:beta_mu[effect[i],i].tolist() for i in effect}
        
        return eff_gamma, eff_mu
    
    def train(self,XX,ytX,XtX,LD,maxite=50,eps=0.01,verbose=False,loss=0.0):
        
        for ite in range(maxite):
            self.infer_q_beta(XX, ytX, XtX, LD)
            ll, mkl, elbo = self.get_elbo()
            if verbose:
                print('*'*70)
                print('Iteration-->{} . Likelihood: {:.2f} . KL: {:.2f} . ELBO: {:.2f}'.format(ite, ll, mkl, elbo))
            if abs(elbo-loss)<eps:
                break
            loss = elbo

parser = argparse.ArgumentParser(description='SparsePro')
parser.add_argument('--ss', type=str, default=None, help='path to summary stats', required=True)
parser.add_argument('--var_Y', type=float, default=None, help='GWAS trait variance', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--LDdir', type=str, default=None, help='path to LD files', required=True)
parser.add_argument('--LDlst', type=str, default=None, help='path to LD list', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', type=bool, default=False, help='options for displaying more information')
parser.add_argument('--tmp', type=bool, default=True, help='options for saving intermediate file')
parser.add_argument('--K', type=int, default=None, help='largest number of effect', required=True)
args = parser.parse_args()

title()

if not os.path.exists(args.save):
    os.makedirs(args.save)

ss = pd.read_csv(args.ss,sep="\s+",dtype={'SNP':str,'BETA':float,'SE':float},index_col=0)
print("summary statistics loaded at {}".format(time.strftime("%Y-%m-%d %H:%M")))

ldlists=pd.read_csv(args.LDlst,header=None)[0]
print("LD list with {} LD blocks loaded\n".format(len(ldlists)))

pip = {}
occ = {}

for i in range(len(ldlists)):
    ld = ldlists[i]
    LD = pd.read_csv(os.path.join(args.LDdir,ld),sep='\t',index_col=0)
    idx = LD.index.intersection(ss.index)
    beta = ss.loc[idx,'BETA'].values
    se = ss.loc[idx,'SE'].values
    XX, XtX, ytX = get_XX_XtX_ytX(LD.values,beta,se,args.var_Y)
    h2_hess=get_HESS_h2_SS(XtX,XX,LD.values,beta,se,args.N,args.var_Y)
    
    print("{} variants loaded from {} with {} variants have matched summary statistics explains {:2.2%} of trait heritability \n".format(LD.shape[1], ld, len(idx), h2_hess))
    model = SparsePro(len(beta),args.K,XX,args.var_Y,h2_hess) 
    model.train(XX, ytX, XtX, LD.values,verbose=args.verbose)
    
    if args.tmp:
        ll,mkl,elbo = model.get_elbo()
        savelist = [h2_hess,model,elbo]
        open_file = open(os.path.join(args.save,'{}.obj'.format(ld)),'wb')
        pickle.dump(savelist,open_file)
        open_file.close()
    
    mcs = model.get_effect_dict()
    # print(otk)
    print('Causal variants for each effect:')
    print({j:[idx[i] for i in mcs[j]] for j in mcs})
    print('Posterior inclusion probabilities for each effect:')
    eff_gamma, eff_mu = model.get_effect_num_dict()
    print(eff_gamma)
    print('Variational effect sizes for each effect:')
    print(eff_mu)
    print()
    
    pip_vec = model.get_PIP().round(4)
    
    for i in range(len(idx)):
        if idx[i] not in pip:
            pip[idx[i]] = pip_vec[i]
            occ[idx[i]] = 1
        else:
            occ[idx[i]] = occ[idx[i]] + 1
            if occ[idx[i]] == 2:
                pip[idx[i]] = pip_vec[i]

allPIP = pd.Series(pip)   
allPIP.to_csv(os.path.join(args.save,"{}.pip".format(args.prefix)),sep='\t',header=None)
print("Statistical fine-mapping finished at {}".format(time.strftime("%Y-%m-%d %H:%M")))
