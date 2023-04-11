
#############################################################################################################################
function isCycle_and_simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j]
  quani=quantum[index]
  qaux[j][1]=qj
  olddx[j][1]=x1j
 
  elapsed = simt - tx[j]
  x[j][0]= xj+elapsed*x1j
  xjaux=x[j][0]
  tx[j]=simt
  u[j][index][1]=u[j][j][1]-aji*qaux[index][1]
  uji=u[j][index][1]
  dxj=aji*qi+ajj*qj+uji
  iscycle=false
  h=0.0;h1=0.0;h2=0.0
  
  if dxj*x1j<0
    qjplus=xjaux+sign(dxj)*quanj
    u[index][j][1]=u[index][index][1]-aij*qj#qaux[j][1]
    uij=u[index][j][1]
    dxi=aii*qi+aij*qjplus+uij
    if dxi*x1i<0
     #=  qaux[j][1]=qj
      olddx[j][1]=x1j =#
     # println("********simul update(val1) double if passed; simt= ",simt)
      iscycle=true  
      println("***************************simul val1 update happened at simt= ",simt)   
      @show index       
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / x1i));h2 = (abs(quanj / x1j));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (0.96*quani / abs(qi - xi));
        h2 = h * (0.96*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 800
          println("maxiter of updateQ      = ",maxIter)
          @show  h, q-x
          
         end
      end

     #=  if 10.0>simt>7.676 && abs(x[2][0]-0.7)<3*quantum[2]
        @show h,h1,h2
        @show xi,qi,quani
        @show xj,qj,quanj
        @show a
          @show u
      end =#

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end
         
function isCycle_and_simulUpdate(::Val{2},index::Int,j::Int,prevStepVal::MVector{T,MVector{B,Float64}},#= direction::MVector{T,Float64}, =# x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}}#= ,olddxSpec::MVector{T,MVector{O,Float64}} =#,tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64,SD::Vector{Vector{Int64}}#= ,qminus::MVector{T,Float64} =##= ,breakloop::MVector{1,Float64}, =##= nextStateTime::MVector{T,Float64} =#)where{T,O,B}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];
    
  e1 = simt - tx[j];e2 = simt - tq[j];e3=simt - tu[j];tu[j]=simt; 
  x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt
   qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  

  xj1=x[j][1]+e1*xj2#= ;olddxSpec[j][1]=xj1 =#;olddx[j][1]=xj1
  # x[j][1]=xj1  # only derx index should be computed using integratestate. other derx should be computed using derx=f()
  

 # ujj=ujj+e1*u[j][j][2]  
 ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  u[j][index][1]=ujj-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uji
  uji=u[j][index][1]

  #u[j][j][2]=xj2-ajj*qj1###################################################-----------------------

  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
  uji2=u[j][index][2] 


  #@show uji2
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  #@show aji,ajj,uji2
  iscycle=false
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
    u[index][j][1]=uii-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2]

    
   
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
  
  
      #|| (dqjplus)*dirj<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
      #β=dxi+sqrt(abs(ddxi)*quani/2)
      #h22=sqrt(2*quani/abs(ddxi))
     # β=dxi+h22*ddxi/2
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) #|| β*direction[index]<0.0

       #=  changexj=abs(dxj-xj1)/(abs(dxj+xj1)/2)+abs(ddxj-xj2)/(abs(ddxj+xj2)/2)
        changexi=abs(dxi-xi1)/(abs(dxi+xi1)/2) + abs(ddxi-xi2)/(abs(ddxi+xi2)/2)
        if changexj>changexi =#
   
        iscycle=true
        isTripleCycle=false

       # flag[index]=1.0
                  #=    for l = 1:length(SD[index])
                      k = SD[index][l] 
                      if k != 0 && k!=index && k!=j && a[index][k]*a[k][index]!=0  
                        akk=a[k][k];aik=a[index][k];aki=a[k][index];xk=x[k][0];qk=q[k][0];qk1=q[k][1];xk1=x[k][1];xk2=2*x[k][2]
                        ukk=u[k][k][1]#;uik=u[index][k][1];uki=u[k][index][1]#;uki2=u[k][index][2]
                        quank=quantum[k];                        
                        e1 = simt - tx[k];e2 = simt - tq[k];e3=simt - tu[k];tx[k]=simt;tu[k]=simt

                        xkaux= x[k](e1);x[k][0]=xkaux
                        qk=qk+e2*qk1  ;qaux[k][1]=qk;q[k][0]=qk;tq[k]=simt
                        

                        newDiff=(xkaux-prevStepVal[k][1])
                        dirk=direction[k]
                        if newDiff*dirk <0.0
                          dirk=-dirk
                        elseif newDiff==0 && dirk!=0.0
                          dirk=0.0
                        elseif newDiff!=0 && dirk==0.0
                          dirk=newDiff
                        else
                        end          
                      
                      
                        xk1=x[k][1]+e1*xk2;olddxSpec[k][1]=xk1;olddx[k][1]=xk1
                        # x[k][1]=xk1  # only derx index should be computed using integratestate. other derx should be computed using derx=f()
                        
                        

                      if abs(qk - xkaux) > 10*quank
                                breakloop[1]+=1.0
                                @show index,k,simt,quank
                                @show qk , xkaux
                                @show x[k][1],xk1,x[k][2],q[k][1]
                                @show nextStateTime[index],nextStateTime[k]
                            
                              end

                        
                        #ukk=ukk+e1*u[k][k][2]  
                        ukk=xk1-akk*qk
                        u[k][k][1]=ukk
                        u[k][index][1]=ukk-aki*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uki
                        uki=u[k][index][1]
                      
                        u[k][k][2]=xk2-akk*qk1###################################################-----------------------
                      
                        u[k][index][2]=u[k][k][2]-aki*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
                        #u[k][index][2]=u[k][k][2]-akk*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
                        uki2=u[k][index][2] 
                      
                      
                        #@show uki2
                        dxk=aki*qi+akk*qaux[k][1]+uki
                        ddxk=aki*qi1+akk*qk1+uki2

                        u[index][k][1]=uii-aik*qaux[k][1]
                        uik=u[index][k][1]
                        u[index][k][2]=u[index][index][2]-aik*qk1#########qaux[k][2] updated in normal Qupdate..ft=20 slightly shifts up
                        uik2=u[index][k][2]



                        qkplus=xkaux-sign(ddxk)*quank
                        h=sqrt(2*quank/abs(ddxk))#2*quantum funny oscillating graph; xk2 vibrating
                        dqkplus=(aki*(qi+h*qi1)+akk*qkplus+uki+h*uki2)/(1-h*akk)
                      
                        dxi=aii*qi+aik*qkplus+uik
                        ddxi=aii*qi1+aik*dqkplus+uik2


                        if (abs(dxk-xk1)>(abs(dxk+xk1)/2) || abs(ddxk-xk2)>(abs(ddxk+xk2)/2))  || (dqkplus)*dirk<0.0 #(dqkplus*qk1)<=0.0 with dir is better since when dir =0 we do not enter
                          #β=dxi+sqrt(abs(ddxi)*quani/2)
                          h3=sqrt(2*quani/abs(ddxi))
                          β=dxi+h3*ddxi/2
                          if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || β*direction[index]<0.0
                                isTripleCycle=true
                                flag[k]=1.0
                                ajk=a[j][k];akj=a[k][j];
                                #@show ajk,akj
                                A=[aii aij aik;aji ajj ajk;aki akj akk]
                                I=[1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0]
                                uijk=uij-aik*qk;ujik=uji-ajk*qk;ukij=uki-akj*qj
                                uijk2=uij2-aik*qk1;ujik2=uji2-ajk*qk1;ukij2=uki2-akj*qj1
                                U=[uijk;ujik;ukij]
                                U2=[uijk2;ujik2;ukij2]
                                X=[xi;xjaux;xkaux]

                               #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj ||abs(qk - xkaux) > 10*quank)
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                                end =#



                                h = ft-simt
                                N=inv(I-h*A)
                                Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))        
                                qi=Q[1]
                                qj=Q[2]
                                qk=Q[3]
                               if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj|| abs(qk - xkaux) > quank) # removing this did nothing...check @btime later
                                  h1 = sqrt(abs(2*quani/ddxi));h2 = sqrt(abs(2*quanj/ddxj));h3 = sqrt(abs(2*quank/ddxk))   #later add derderX =1e-12 when x2==0
                                  h=min(h1,h2,h3)
                                  N=inv(I-h*A)
                                  Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
                                  qi=Q[1]
                                  qj=Q[2]
                                  qk=Q[3]
                                end
                               
                               #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj ||abs(qk - xkaux) > 10*quank)
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                                end =#

                                maxIter=600
                                while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj ||abs(qk - xkaux) > quank) && (maxIter>0)
                                  maxIter-=1
                                  h1 = h * (0.99*quani / abs(qi - xi));
                                  h2 = h * (0.99*quanj / abs(qj - xjaux));
                                  h3 = h * (0.99*quank / abs(qk - xkaux));
                                  h=min(h1,h2,h3)
                                  N=inv(I-h*A)
                                  Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+( X+(h*U+(h*h*U2/2))) )
                                  qi=Q[1]
                                  qj=Q[2]
                                  qk=Q[3]
                                end
                                Q1=N*(A*Q+U+h*U2)
                                q[index][1]=Q1[1]# store back helper vars
                                q[j][1]=Q1[2]
                                q[k][1]=Q1[3]
                              #  @show simt
                                  ######################
                                  if maxIter < 20
                                    @show maxIter,h
                                    @show simt
                                    @show a
                                 end
                                 #= if h <0.0000001
                                  @show maxIter,index,j,k
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                               end =#
                                 # tq[j]=simt
                                 # tq[k]=simt
                          end #end inner if check triple
                        end #end outer if check triple
                      end# end  if k != 0
                    end#end for to find k 
 =#
       

          if !isTripleCycle
            h = ft-simt
        #= A=[aii aij;aji ajj]
        I=[1.0 0.0;0.0 1.0]
        U=[uij;uji]
        U2=[uij2;uji2]
        X=[xi;xjaux] =#
        Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       # N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
       αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
       αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
       αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
     #  NA=[αii αij;αji αjj]
       βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
       βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
       βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
       βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
    #   M=[βii βij;βji βjj]
       Δ2=βii*βjj-βij*βji
      # invM=[βjj -βij;-βji βii]/Δ2
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]
      #  upart2=[uij+h*uij2;uji+h*uji2]
       # apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]
        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
     #   apart1N=[λii λij;λji λjj]/Δ1
       # part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
      #  part=[parti;partj]
        # part=[part1[1]+xpart2[1];part1[2]+xpart2[2]]
       #=  qi=(βjj*parti-βij*partj)/Δ2
        qj=(βii*partj-βji*parti)/Δ2 =#
        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
        # Q=invM*(part1+(xpart2))
       
       # Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(xpart2)) 
        
      
       #=  qi=Q[1]
        qj=Q[2] =#
       #= Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
       αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
       αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
       αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
       βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
       βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
       βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
       βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
       Δ2=βii*βjj-βij*βji
       γi=(h*h*aii/2-h)*(uij+h*uij2)+(uji+h*uji2)*h*h*aij/2+(xi+h*uij+h*h*uij2/2)
       γj=h*h*aji*(uij+h*uij2)/2+(h*h*ajj/2-h)*(uji+h*uji2)+(xjaux+h*uji+h*h*uji2/2)
       qi=(βjj*γi-βij*γj)/Δ2
       qj=(βii*γj-βji*γi)/Δ2 =#




        if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) # removing this did nothing...check @btime later
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          #N=inv(I-h*A)
          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
    #    N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
        αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
        αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
        αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
       # NA=[αii αij;αji αjj]
        βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
        βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
        βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
        βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
     #   M=[βii βij;βji βjj]
         #N=inv(I-h*A)
         Δ2=βii*βjj-βij*βji
     #  invM=[βjj -βij;-βji βii]/Δ2
        #N=inv(I-h*A)
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]

      #  upart2=[uij+h*uij2;uji+h*uji2]
      #  apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]


        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
       # apart1N=[λii λij;λji λjj]/Δ1
      #  part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
       # part=[parti;partj]
        # part=[part1[1]+xpart2[1];part1[2]+xpart2[2]]
         # Q=invM*(part)
        # invM2=[βjj -βij;-βji βii]/Δ2
        
         qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
        #Q=invM*(part1+(xpart2))
       # Q=invM*(apart1N*(upart2)+(xpart2))
       # Q=invM*((apart1)*N*(upart2)+(xpart2))

        #Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(xpart2)) 
        #Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))   
        # Q=inv(M)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2)))) 
        #  Q=inv(I-h*A+h*NA-h*h*A*NA/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
         # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)  
         #=  qi=Q[1]
          qj=Q[2] =#
         #=  Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
          αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
          αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
          αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
          βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
          βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
          βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
          βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
          Δ2=βii*βjj-βij*βji
          γi=(h*h*aii/2-h)*(uij+h*uij2)+(uji+h*uji2)*h*h*aij/2+(xi+h*uij+h*h*uij2/2)
          γj=h*h*aji*(uij+h*uij2)/2+(h*h*ajj/2-h)*(uji+h*uji2)+(xjaux+h*uji+h*h*uji2/2)
          qi=(βjj*γi-βij*γj)/Δ2
          qj=(βii*γj-βji*γi)/Δ2 =#
        end
       
        maxIter=600
        while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (0.98*quani / abs(qi - xi));
          h2 = h * (0.98*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          #N=inv(I-h*A)
          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        #  N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
          αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
          αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
          αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
          αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
         # NA=[αii αij;αji αjj]
          βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
          βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
          βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
          βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
        #  M=[βii βij;βji βjj]
           #N=inv(I-h*A)
           Δ2=βii*βjj-βij*βji
     #  invM=[βjj -βij;-βji βii]/Δ2
        #N=inv(I-h*A)
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]

      #  upart2=[uij+h*uij2;uji+h*uji2]
       # apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]


        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
      #  apart1N=[λii λij;λji λjj]/Δ1
       # part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
      
        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)


       
          
        end
    
        if maxIter < 20
           @show maxIter
           @show simt
           @show a
        end
       if debug @show maxIter  end

     
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
         q1parti=aii*qi+aij*qj+uij+h*uij2
         q1partj=aji*qi+ajj*qj+uji+h*uji2
      
         q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
        q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
           
          end #end double cycle done
       # end#if changexj>changexi
      end #end second dependecy check
    end # end outer dependency check
  return iscycle
end
function nisCycle_and_simulUpdate(::Val{2},index::Int,j::Int,prevStepVal::MVector{T,MVector{B,Float64}},direction::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},olddxSpec::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64,SD::Vector{Vector{Int64}},qminus::MVector{T,Float64}#= ,breakloop::MVector{1,Float64}, =##= nextStateTime::MVector{T,Float64} =#)where{T,O,B}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];
    
  e1 = simt - tx[j];e2 = simt - tq[j];e3=simt - tu[j];tu[j]=simt; 
  x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt
  qminus[j]=qj
  qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  

  newDiff=(xjaux-prevStepVal[j][1])
  dirj=direction[j]
  if newDiff*dirj <0.0
    dirj=-dirj
 
  elseif newDiff==0 && dirj!=0.0
    dirj=0.0
   
  elseif newDiff!=0 && dirj==0.0
    dirj=newDiff
  else
 
  end          
  

  xj1=x[j][1]+e1*xj2;olddxSpec[j][1]=xj1;olddx[j][1]=xj1
  # x[j][1]=xj1  # only derx index should be computed using integratestate. other derx should be computed using derx=f()
  

  #ujj=ujj+e1*u[j][j][2]  
  ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  u[j][index][1]=ujj-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uji
  uji=u[j][index][1]

  u[j][j][2]=xj2-ajj*qj1###################################################-----------------------

  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
  uji2=u[j][index][2] 


  #@show uji2
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  #@show aji,ajj,uji2
  iscycle=false
  
    u[index][j][1]=uii-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2]

    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
   
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
  
  
    if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || (dqjplus)*dirj<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
      #β=dxi+sqrt(abs(ddxi)*quani/2)
      h2=sqrt(2*quani/abs(ddxi))
      β=dxi+h2*ddxi/2
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || β*direction[index]<0.0

       #=  changexj=abs(dxj-xj1)/(abs(dxj+xj1)/2)+abs(ddxj-xj2)/(abs(ddxj+xj2)/2)
        changexi=abs(dxi-xi1)/(abs(dxi+xi1)/2) + abs(ddxi-xi2)/(abs(ddxi+xi2)/2)
        if changexj>changexi =#
   
        iscycle=true
        isTripleCycle=false

       # flag[index]=1.0
                  #=    for l = 1:length(SD[index])
                      k = SD[index][l] 
                      if k != 0 && k!=index && k!=j && a[index][k]*a[k][index]!=0  
                        akk=a[k][k];aik=a[index][k];aki=a[k][index];xk=x[k][0];qk=q[k][0];qk1=q[k][1];xk1=x[k][1];xk2=2*x[k][2]
                        ukk=u[k][k][1]#;uik=u[index][k][1];uki=u[k][index][1]#;uki2=u[k][index][2]
                        quank=quantum[k];                        
                        e1 = simt - tx[k];e2 = simt - tq[k];e3=simt - tu[k];tx[k]=simt;tu[k]=simt

                        xkaux= x[k](e1);x[k][0]=xkaux
                        qk=qk+e2*qk1  ;qaux[k][1]=qk;q[k][0]=qk;tq[k]=simt
                        

                        newDiff=(xkaux-prevStepVal[k][1])
                        dirk=direction[k]
                        if newDiff*dirk <0.0
                          dirk=-dirk
                        elseif newDiff==0 && dirk!=0.0
                          dirk=0.0
                        elseif newDiff!=0 && dirk==0.0
                          dirk=newDiff
                        else
                        end          
                      
                      
                        xk1=x[k][1]+e1*xk2;olddxSpec[k][1]=xk1;olddx[k][1]=xk1
                        # x[k][1]=xk1  # only derx index should be computed using integratestate. other derx should be computed using derx=f()
                        
                        

                      if abs(qk - xkaux) > 10*quank
                                breakloop[1]+=1.0
                                @show index,k,simt,quank
                                @show qk , xkaux
                                @show x[k][1],xk1,x[k][2],q[k][1]
                                @show nextStateTime[index],nextStateTime[k]
                            
                              end

                        
                        #ukk=ukk+e1*u[k][k][2]  
                        ukk=xk1-akk*qk
                        u[k][k][1]=ukk
                        u[k][index][1]=ukk-aki*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uki
                        uki=u[k][index][1]
                      
                        u[k][k][2]=xk2-akk*qk1###################################################-----------------------
                      
                        u[k][index][2]=u[k][k][2]-aki*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
                        #u[k][index][2]=u[k][k][2]-akk*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
                        uki2=u[k][index][2] 
                      
                      
                        #@show uki2
                        dxk=aki*qi+akk*qaux[k][1]+uki
                        ddxk=aki*qi1+akk*qk1+uki2

                        u[index][k][1]=uii-aik*qaux[k][1]
                        uik=u[index][k][1]
                        u[index][k][2]=u[index][index][2]-aik*qk1#########qaux[k][2] updated in normal Qupdate..ft=20 slightly shifts up
                        uik2=u[index][k][2]



                        qkplus=xkaux-sign(ddxk)*quank
                        h=sqrt(2*quank/abs(ddxk))#2*quantum funny oscillating graph; xk2 vibrating
                        dqkplus=(aki*(qi+h*qi1)+akk*qkplus+uki+h*uki2)/(1-h*akk)
                      
                        dxi=aii*qi+aik*qkplus+uik
                        ddxi=aii*qi1+aik*dqkplus+uik2


                        if (abs(dxk-xk1)>(abs(dxk+xk1)/2) || abs(ddxk-xk2)>(abs(ddxk+xk2)/2))  || (dqkplus)*dirk<0.0 #(dqkplus*qk1)<=0.0 with dir is better since when dir =0 we do not enter
                          #β=dxi+sqrt(abs(ddxi)*quani/2)
                          h3=sqrt(2*quani/abs(ddxi))
                          β=dxi+h3*ddxi/2
                          if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || β*direction[index]<0.0
                                isTripleCycle=true
                                flag[k]=1.0
                                ajk=a[j][k];akj=a[k][j];
                                #@show ajk,akj
                                A=[aii aij aik;aji ajj ajk;aki akj akk]
                                I=[1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0]
                                uijk=uij-aik*qk;ujik=uji-ajk*qk;ukij=uki-akj*qj
                                uijk2=uij2-aik*qk1;ujik2=uji2-ajk*qk1;ukij2=uki2-akj*qj1
                                U=[uijk;ujik;ukij]
                                U2=[uijk2;ujik2;ukij2]
                                X=[xi;xjaux;xkaux]

                               #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj ||abs(qk - xkaux) > 10*quank)
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                                end =#



                                h = ft-simt
                                N=inv(I-h*A)
                                Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))        
                                qi=Q[1]
                                qj=Q[2]
                                qk=Q[3]
                               if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj|| abs(qk - xkaux) > quank) # removing this did nothing...check @btime later
                                  h1 = sqrt(abs(2*quani/ddxi));h2 = sqrt(abs(2*quanj/ddxj));h3 = sqrt(abs(2*quank/ddxk))   #later add derderX =1e-12 when x2==0
                                  h=min(h1,h2,h3)
                                  N=inv(I-h*A)
                                  Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
                                  qi=Q[1]
                                  qj=Q[2]
                                  qk=Q[3]
                                end
                               
                               #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj ||abs(qk - xkaux) > 10*quank)
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                                end =#

                                maxIter=600
                                while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj ||abs(qk - xkaux) > quank) && (maxIter>0)
                                  maxIter-=1
                                  h1 = h * (0.99*quani / abs(qi - xi));
                                  h2 = h * (0.99*quanj / abs(qj - xjaux));
                                  h3 = h * (0.99*quank / abs(qk - xkaux));
                                  h=min(h1,h2,h3)
                                  N=inv(I-h*A)
                                  Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+( X+(h*U+(h*h*U2/2))) )
                                  qi=Q[1]
                                  qj=Q[2]
                                  qk=Q[3]
                                end
                                Q1=N*(A*Q+U+h*U2)
                                q[index][1]=Q1[1]# store back helper vars
                                q[j][1]=Q1[2]
                                q[k][1]=Q1[3]
                              #  @show simt
                                  ######################
                                  if maxIter < 20
                                    @show maxIter,h
                                    @show simt
                                    @show a
                                 end
                                 #= if h <0.0000001
                                  @show maxIter,index,j,k
                                  @show simt,h
                                  @show abs(qi - xi)
                                  @show abs(qj - xjaux)
                                  @show abs(qk - xkaux)
                               end =#
                                 # tq[j]=simt
                                 # tq[k]=simt
                          end #end inner if check triple
                        end #end outer if check triple
                      end# end  if k != 0
                    end#end for to find k 
 =#
       

          if !isTripleCycle
            h = ft-simt
        #= A=[aii aij;aji ajj]
        I=[1.0 0.0;0.0 1.0]
        U=[uij;uji]
        U2=[uij2;uji2]
        X=[xi;xjaux] =#
        Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       # N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
       αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
       αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
       αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
     #  NA=[αii αij;αji αjj]
       βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
       βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
       βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
       βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
    #   M=[βii βij;βji βjj]
       Δ2=βii*βjj-βij*βji
      # invM=[βjj -βij;-βji βii]/Δ2
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]
      #  upart2=[uij+h*uij2;uji+h*uji2]
       # apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]
        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
     #   apart1N=[λii λij;λji λjj]/Δ1
       # part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
      #  part=[parti;partj]
        # part=[part1[1]+xpart2[1];part1[2]+xpart2[2]]
       #=  qi=(βjj*parti-βij*partj)/Δ2
        qj=(βii*partj-βji*parti)/Δ2 =#
        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
        # Q=invM*(part1+(xpart2))
       
       # Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(xpart2)) 
        
      
       #=  qi=Q[1]
        qj=Q[2] =#
       #= Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
       αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
       αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
       αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
       αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
       βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
       βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
       βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
       βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
       Δ2=βii*βjj-βij*βji
       γi=(h*h*aii/2-h)*(uij+h*uij2)+(uji+h*uji2)*h*h*aij/2+(xi+h*uij+h*h*uij2/2)
       γj=h*h*aji*(uij+h*uij2)/2+(h*h*ajj/2-h)*(uji+h*uji2)+(xjaux+h*uji+h*h*uji2/2)
       qi=(βjj*γi-βij*γj)/Δ2
       qj=(βii*γj-βji*γi)/Δ2 =#




        if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) # removing this did nothing...check @btime later
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          #N=inv(I-h*A)
          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
    #    N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
        αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
        αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
        αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
       # NA=[αii αij;αji αjj]
        βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
        βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
        βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
        βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
     #   M=[βii βij;βji βjj]
         #N=inv(I-h*A)
         Δ2=βii*βjj-βij*βji
     #  invM=[βjj -βij;-βji βii]/Δ2
        #N=inv(I-h*A)
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]

      #  upart2=[uij+h*uij2;uji+h*uji2]
      #  apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]


        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
       # apart1N=[λii λij;λji λjj]/Δ1
      #  part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
       # part=[parti;partj]
        # part=[part1[1]+xpart2[1];part1[2]+xpart2[2]]
         # Q=invM*(part)
        # invM2=[βjj -βij;-βji βii]/Δ2
        
         qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
        #Q=invM*(part1+(xpart2))
       # Q=invM*(apart1N*(upart2)+(xpart2))
       # Q=invM*((apart1)*N*(upart2)+(xpart2))

        #Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(xpart2)) 
        #Q=invM*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))   
        # Q=inv(M)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2)))) 
        #  Q=inv(I-h*A+h*NA-h*h*A*NA/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
         # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)  
         #=  qi=Q[1]
          qj=Q[2] =#
         #=  Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
          αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
          αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
          αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
          βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
          βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
          βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
          βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
          Δ2=βii*βjj-βij*βji
          γi=(h*h*aii/2-h)*(uij+h*uij2)+(uji+h*uji2)*h*h*aij/2+(xi+h*uij+h*h*uij2/2)
          γj=h*h*aji*(uij+h*uij2)/2+(h*h*ajj/2-h)*(uji+h*uji2)+(xjaux+h*uji+h*h*uji2/2)
          qi=(βjj*γi-βij*γj)/Δ2
          qj=(βii*γj-βji*γi)/Δ2 =#
        end
       
        maxIter=600
        while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (0.98*quani / abs(qi - xi));
          h2 = h * (0.98*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          #N=inv(I-h*A)
          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        #  N=[1-h*ajj h*aij;h*aji 1-h*aii]/Δ1
          αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
          αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
          αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
          αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
         # NA=[αii αij;αji αjj]
          βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
          βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
          βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
          βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
        #  M=[βii βij;βji βjj]
           #N=inv(I-h*A)
           Δ2=βii*βjj-βij*βji
     #  invM=[βjj -βij;-βji βii]/Δ2
        #N=inv(I-h*A)
      #  xpart2=[xi+h*uij+h*h*uij2/2;xjaux+h*uji+h*h*uji2/2]

      #  upart2=[uij+h*uij2;uji+h*uji2]
       # apart1=[h*h*aii/2-h h*h*aij/2;h*h*aji/2 h*h*ajj/2-h]


        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
      #  apart1N=[λii λij;λji λjj]/Δ1
       # part1=[λii*(uij+h*uij2)+λij*(uji+h*uji2);λji*(uij+h*uij2)+λjj*(uji+h*uji2)]/Δ1
        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
      
        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)


       
          
        end
    
        if maxIter < 20
           @show maxIter
           @show simt
           @show a
        end
       if debug @show maxIter  end

     
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
         q1parti=aii*qi+aij*qj+uij+h*uij2
         q1partj=aji*qi+ajj*qj+uji+h*uji2
      
         q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
        q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
           
          end #end double cycle done
       # end#if changexj>changexi
      end #end second dependecy check
    end # end outer dependency check
  return iscycle
end

function isCycle_and_simulUpdate(::Val{3},index::Int,j::Int,prevStepVal::MVector{T,MVector{B,Float64}},direction::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},olddxSpec::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64,SD::Vector{Vector{Int64}}#= ,flag::MVector{T,Float64} =#)where{T,O,B}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qi2=2*q[index][2];qj1=q[j][1];qj2=2*q[j][2]
  uii=u[index][index][1];uij=u[index][j][1];ujj=u[j][j][1]#;uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];xi1=x[index][1];xi2=2*x[index][2];xi3=6*x[index][3];xj1=x[j][1];xj2=2*x[j][2];xj3=6*x[j][3]
  e1 = simt - tx[j];e2 = simt - tq[j];e3=simt - tu[j]; tx[j]=simt;tu[j]=simt
  x[j][0]= x[j](e1);xjaux=x[j][0]


  xj1=xj1+e1*xj2+e1*e1*xj3/2;olddxSpec[j][1]=xj1;olddx[j][1]=xj1; xj2=xj2+e1*xj3
  qj=qj+e2*qj1+e2*e2*qj2/2  ;qj1=qj1+e2*qj2; qaux[j][1]=qj ;qaux[j][2]=qj1

  newDiff=(x[j][0]-prevStepVal[j][1])
  dir=direction[j]
  if newDiff*dir <0.0
    direction[j]=-dir
   # @show simt
    #@show k,direction[k]
  elseif newDiff==0 && dir!=0.0
    direction[j]=0.0
    #@show simt
    #@show k,direction[k]
  elseif newDiff!=0 && dir==0.0
    direction[j]=newDiff
  else
  #do not update direction
  end          




  #ujj=ujj+e1*u[j][j][2]+e1*e1*u[j][j][3]/2  
  ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  #u[j][j][2]=u[j][j][2]+e1*u[j][j][3]
  u[j][j][2]=xj2-ajj*qj1
  u[j][j][3]=xj3-ajj*qj2###################################
  u[j][index][1]=ujj-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uji
  uji=u[j][index][1]
  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
  uji2=u[j][index][2]
   u[j][index][3]=u[j][j][3]-aji*qaux[index][3]#
   uji3=u[j][index][3]

  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  dddxj=aji*qi2+ajj*qj2+uji3
  #@show aji,ajj,uji,uji2,uji3
  iscycle=false

    u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1
    uij2=u[index][j][2]
    u[index][j][3]=u[index][index][3]-aij*qj2
    uij3=u[index][j][3]



 # if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2) || abs(dddxj-xj3)>(abs(dddxj+xj3)/2))  
    
#=     if abs(dxj-xj1)>(abs(dxj+xj1)/2)
      println("1st der is culprit")
      @show dxj,xj1
    end
    if abs(ddxj-xj2)>(abs(ddxj+xj2)/2)
      println("2nd der is culprit")
      @show ddxj,xj2
    end
    if abs(dddxj-xj3)>(abs(dddxj+xj3)/2)
      println("3rd der is culprit")
      @show dddxj,xj3
    end =#



    h1=cbrt(abs(6*quanj/dddxj))
    qjplus=xjaux+h1*h1*h1*dddxj/6
    #@show h,qjplus,xjaux
    λ=ajj*qjplus+aji*qi+uji+h1*(aji*qi1+uji2)+h1*h1*(aji*qi2+uji3)/2
    α=(h1-h1*h1/2)*(aji*(qi1+h1*qi2)+uji2+h1*uji3)/(1-h1*ajj)
    β=1-h1*ajj+(h1-h1*h1/2)*ajj/(1-h1*ajj)
    dqjplus=(λ-α)/β
    ddqjplus=(aji*(qi1+h1*qi2)+ajj*dqjplus+uji2+h1*uji3)/(1-h1*ajj)
   ### 
    #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
   
   #=  u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2] =#
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
    dddxi=aii*qi2+aij*ddqjplus+uij3
   #=  if debug
    println("----in case the first if passed isCycle function modified these: ")
    @show uij
    @show uij2
    @show dxi,xi1
    @show ddxi,xi2
    end =#
 #=  #dqjplus-qj1 is enough: normally (qjplus+h*dqjlpus-qj+h*qj1) . dqjplus-qj1 is better than dqjplus alone cuz the case dqjplus=0
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || (dqjplus+h1*ddqjplus)*direction[j]<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    h2=cbrt(abs(6*quani/dddxi))
    β=dxi+ddxi*h2/2+dddxi*h2*h2/6
    if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || β*direction[index]<0.0 =#

      if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2) || abs(dddxj-xj3)>(abs(dddxj+xj3)/2))  || (dqjplus+h1*ddqjplus)*direction[j]<0.0#(dqjplus+h1*ddqjplus-qj1-e1*qj2)*direction[j]<0.0
        h2=cbrt(abs(6*quani/dddxi))
      h2=cbrt(abs(6*quani/dddxi))
      β=dxi+ddxi*h2/2+dddxi*h2*h2/6
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)  || abs(dddxi-xi3)>(abs(dddxi+xi3)/2))|| β*direction[index]<0.0

   # if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)  || abs(dddxi-xi3)>(abs(dddxi+xi3)/2))

        iscycle=true
  

        A=[aii aij;aji ajj]
        I=[1.0 0.0;0.0 1.0]
        U=[uij;uji]
        U2=[uij2;uji2]
        U3=[uij3;uji3]
        X=[xi;xjaux]


                                       h = ft-simt
                                        N=inv(I-h*A)
                                        R=h*I-(h*h/2)*A+(((h*h/2)*I-(h*h*h/6)*A)*N*A)
                                        L=inv(I-h*A+h*N*A-(h*h/2)*A*N*A)
                                        M=((h*h/2)*A-h*I)*N*(U2+h*U3)
                                        P=X+(h*U+(h*h/2)*U2+(h*h*h/6)*U3)
                                        S=-R*L*(M+U+h*U2+(h*h/2)*U3)-(((h*h/2)*I-(h*h*h/6)*A)*N*(U2+h*U3))+P
                                        Q=inv(I-h*A+R*L*A)*S
                                        
                                        qi=Q[1]
                                        qj=Q[2]
                                      
                                        if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) # removing this did nothing...check @btime later
                                          h1 = cbrt(abs(6*quani/dddxi));h2 = cbrt(abs(6*quanj/dddxj));   #later add derderX =1e-12 when x2==0
                                          h=min(h1,h2)
                                          N=inv(I-h*A)
                                          R=h*I-(h*h/2)*A+(((h*h/2)*I-(h*h*h/6)*A)*N*A)
                                        L=inv(I-h*A+h*N*A-(h*h/2)*A*N*A)
                                        M=((h*h/2)*A-h*I)*N*(U2+h*U3)
                                        P=X+(h*U+(h*h/2)*U2+(h*h*h/6)*U3)
                                        S=-R*L*(M+U+h*U2+(h*h/2)*U3)-(((h*h/2)*I-(h*h*h/6)*A)*N*(U2+h*U3))+P
                                        Q=inv(I-h*A+R*L*A)*S
                                        # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)  
                                          qi=Q[1]
                                          qj=Q[2]
                                        end
                                       #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj) 
                                          @show xi,qi
                                          @show xjaux,qj
                                          @show h
                                    end =#
                                        maxIter=600000
                                        while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0) && (h!=0.0)
                                          maxIter-=1
                                          h1 = h * (0.98*quani / abs(qi - xi));
                                          h2 = h * (0.98*quanj / abs(qj - xjaux));
                                          h=min(h1,h2)
                                        # h=h-0.00001
                                          N=inv(I-h*A)
                                          R=h*I-(h*h/2)*A+(((h*h/2)*I-(h*h*h/6)*A)*N*A)
                                          L=inv(I-h*A+h*N*A-(h*h/2)*A*N*A)
                                          M=((h*h/2)*A-h*I)*N*(U2+h*U3)
                                          P=X+(h*U+(h*h/2)*U2+(h*h*h/6)*U3)
                                          S=-R*L*(M+U+h*U2+(h*h/2)*U3)-(((h*h/2)*I-(h*h*h/6)*A)*N*(U2+h*U3))+P
                                          Q=inv(I-h*A+R*L*A)*S
                                          
                                          qi=Q[1]
                                          qj=Q[2]
                                          
                                        end
                                       #=  if (abs(qi - xi) > 10*quani || abs(qj - xjaux) > 10*quanj) 
                                        println("after while")
                                          @show xi,qi
                                          @show xjaux,qj
                                          @show h
                                    end =#
                                      
                                      
                                                                                       if maxIter < 20
                                                                                        println("maxtiter simult_val{3}")
                                                                                        @show maxIter
                                                                                        
                                                                                      end
                                                                          #=if debug @show maxIter  end =#
                                        q[index][0]=qi# store back helper vars
                                        q[j][0]=qj
                                      
                                      
                                        Q1=L*(M+A*Q+U+h*U2+(h*h/2)*U3)
                                        q[index][1]=Q1[1]# store back helper vars
                                        q[j][1]=Q1[2]
       #=  h = ft-simt
        N=inv(I-h*A)
        Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))    
       # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)      
        qi=Q[1]
        qj=Q[2]
       


        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) # removing this did nothing...check @btime later
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
         # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)  
          qi=Q[1]
          qj=Q[2]
        end
        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
          @show xi,qi
          @show xjaux,qj
          @show h
    end
        maxIter=600
        while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (1.96*quani / abs(qi - xi));
          h2 = h * (1.96*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          N=inv(I-h*A)
         
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+( X+(h*U+(h*h*U2/2))) )
           
          
          qi=Q[1]
          qj=Q[2]
          
        end
  
  
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        Q1=N*(A*Q+U+h*U2)
        q[index][1]=Q1[1]# store back helper vars
        q[j][1]=Q1[2]
 =#
        Q2=N*(A*Q1+U2+h*U3)

        q[index][2]=Q2[1]/2# store back helper vars: /2 for taylor standard storage
        q[j][2]=Q2[2]/2

        #= if  abs(qi+h*Q1[1]-xi)<1e-8  || abs(qj+h*Q1[2]-xjaux)<1e-8
          println("next x is about the same this x")
          @show h
          @show qi+h*Q1[1],xi
          @show qj+h*Q1[2],xjaux
          @show simt
          @show j
        end
        if abs(qi-xi)<1e-12 || abs(qj-xjaux)<1e-12
          println("could not throw q")
          @show qi,xi
          @show qj,xjaux
          @show simt
          @show j
        end
        if debug
          println("h used was = ",h) 
          @show q[index]
          @show q[j]       
          println("end of iscycle function")
         end =#
        tq[j]=simt
      end #end second dependecy check
  end # end outer dependency check
  #= if debug  
    println("end of iscycle function")
   end =#
  return iscycle
end
#######################################################################################################################################################
function Liqss_simulreComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
  dt=0.0
  q=qv[i][0]
  x=xv[i][0]
  if xv[i][1] !=0.0  #&& abs(q-x)>quantum[i]/10
      if abs(q-x)>quantum[i]/10 # if q and x are very close then equilibrium no need to enter to compute next
          dt=(q-x)/xv[i][1]
          if dt>0.0
              nextTime[i]=currentTime+dt# later guard against very small dt
          elseif dt<0.0
              if xv[i][1]>0.0  
                  nextTime[i]=currentTime+(quantum[i])/xv[i][1]
              else
                  nextTime[i]=currentTime+(-quantum[i])/xv[i][1]
              end
             #=  if nextTime[i] < currentTime 
                  println("liqss1 simulrecompute: x heading neither towards nor 2delta away! at time= ",currentTime) 
                  nextTime[i]=currentTime+1e-1
              end =#
          else#q=x equilibruim point reached#not needed with if abs(q-x)>quantum[i]/10
              nextTime[i]=Inf
              println("simulrecompval1 dt==0")
          end
          #= if 10.0>currentTime>7.676
              println(" nextTime[$i]= ", nextTime[i])
            end =#
      else
          nextTime[i]=Inf
          if nextTime[1]==Inf && nextTime[2]==Inf
          println("!!!!!!!!!!!!!!!!!!!!!!!!liqss1 simulrecompute: abs(q-x)<quantum[$i]/10!   equilibrium at time? ",currentTime) 
          end
      end
  else
      nextTime[i]=Inf
  end
end



function Liqss_simulreComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
  q=qv[i][0]
  x=xv[i][0]
  q1=qv[i][1]
  x1=xv[i][1]
  x2=xv[i][2]
  if abs(q-x)>quantum[i]/10
  coef=@SVector [q-x, q1-x1,-x2]#
      nextTime[i]=currentTime + minPosRoot(coef, Val(2))
      if q-x >0
          coef=setindex(coef, -quantum[i],1)
          timetemp = currentTime + minPosRoot(coef, Val(2))
          if timetemp < nextTime[i] 
              nextTime[i]=timetemp
          end
      elseif  q-x <0
          coef=setindex(coef, quantum[i],1)
          timetemp = currentTime + minPosRoot(coef, Val(2))
          if timetemp < nextTime[i] 
              nextTime[i]=timetemp
          end
      else
          nextTime[i]=Inf #
          println("equilibrium reached")
      end
    else
      nextTime[i]=Inf
      if nextTime[1]==Inf && nextTime[2]==Inf
      println("!!!!!!!!!!!!!!!!!!!!!!!!liqss2 simulrecompute: abs(q-x)<quantum[$i]/10!   equilibrium at time? ",currentTime) 
      end
  end
  #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
      time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
      nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
  end  =#   
end


#######################################################################################################################################################
function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::Vector{Vector{Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
      diffQ=(q[index][0]-qaux[index][1])
  
        if diffQ!=0
        a[j][index]=(x[j][1]-olddx[j][1])/diffQ
        else
        a[j][index]=0.0
        end
     
        return nothing
end
function updateOtherApprox(::Val{2},k::Int,j::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
  diffQ=q[j][0]-qaux[j][1]
  # println("aji before updateoher= ",a[j][index])
  if diffQ != 0.0
      a[k][j]=(x[k][1]-olddx[k][1])/diffQ
  else
      a[k][j]=0.0
  end
 
     u[k][j][1]=x[k][1]-a[k][k]*q[k][0]-a[k][j]*q[j][0]
   u[k][j][2]=2*x[k][2]-a[k][k]*q[k][1]-a[k][j]*q[j][1]
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
  return nothing
end

function updateOtherApprox(::Val{3},k::Int,j::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
  diffQ=q[j][0]-qaux[j][1]
  # println("aji before updateoher= ",a[j][index])
  if diffQ != 0.0
      a[k][j]=(x[k][1]-olddx[k][1])/diffQ
  else
      a[k][j]=0.0
  end
 
     u[k][j][1]=x[k][1]-a[k][k]*q[k][0]-a[k][j]*q[j][0]
   u[k][j][2]=2*x[k][2]-a[k][k]*q[k][1]-a[k][j]*q[j][1]
   u[k][j][3]=6*x[k][3]-a[k][k]*2*q[k][2]-a[k][j]*2*q[j][2]
  #  println("u inside updateOther after update= ",u[j][index])
  #tu[k]=simt
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
  return nothing
end
