 
if k in jac[k]
    qctemp=q[k][0];q[k][0]=qminus[k]# i want only the effect of qc on acc: remove influence of index and j
    clearCache(taylorOpsCache,cacheSize);f(k,-1,-1,q,d,t,taylorOpsCache);
    olddx[k][1]=taylorOpsCache[1][0]# acc =(dxc-oldxc)/(qc-qcminus)
    q[k][0]=qctemp
    nupdateLinearApprox(Val(O),k,x,q,a,u,qminus[k],olddx,tu,simt)# acc =(dxc-oldxc)/(qc-qcminus)
  else
    nupdateUaNull(Val(O),k,x,u,tu,simt)# acc==0
  end