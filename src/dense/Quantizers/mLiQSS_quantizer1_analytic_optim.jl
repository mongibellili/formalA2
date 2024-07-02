function misCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,rejected,simpleCase,temp2,temp3)
  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index,simt)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,d,cacheA,j,j,simt)
  ajj=cacheA[1]
  α=-aii-ajj
    β=aii*ajj-aij*aji
    if α<0 || β<0
      rejected[1]+=1
      return false
    end

   xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  #qaux[j][1]=qj;
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;

  xj=x[j][0]
  tx[j]=simt

 qiminus=qaux[index][1]
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qiminus
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qiminus-aij*qj
    iscycle=false
    dxj=aji*qi+ajj*qj+uji #only future qi   #emulate fj
   
    dxithrow=aii*qi+aij*qj+uij #only future qi
                                                      
  qjplus=xj+sign(dxj)*quanj  #emulate updateQ(j)...

    dxi=aii*qi+aij*qjplus+uij #both future qi & qj   #emulate fi
  
  ########old condition:Union 
     #= if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
    if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
        iscycle=true
    end
  end   =#

    ########condition:Union i union
    if (abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0)
      cancelCriteria=1e-6*quani
      if abs(dxi)>3*abs(dxithrow) || abs(dxi)*3<abs(dxithrow) ||  (dxi*dxithrow)<0.0  
          iscycle=true
          if abs(dxi)<cancelCriteria && abs(dxithrow)<cancelCriteria #significant change is relative to der values. cancel if all values are small
            iscycle=false
           # println("cancel dxi & dxthr small")
          end
      end
      if  abs(dxi)>10*abs(ẋi) || abs(dxi)*10<abs(ẋi) #||  (dxi*ẋi)<0.0 
        iscycle=true
        if abs(dxi)<cancelCriteria && abs(ẋi)<cancelCriteria #significant change is relative to der values. cancel if all values are small
          iscycle=false
         # println("cancel dxi & ẋi small")
        end
      end
      
      if abs(dxj)<1e-6*quanj && abs(ẋj)<1e-6*quanj #significant change is relative to der values. cancel if all values are small
        iscycle=false
        #println("cancel dxj & ẋj small")
      end

    end   
   

   # iscycle=false
  if iscycle

 
   
      
       #find positive zeros f=+-Δ
      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij
    

      if bi*ci>0.0 && bj*cj>0.0
        simpleCase[1]+=1
        hi=bcPOS(bi,ci,quani,βi,αi)
        hj=bcPOS(bj,cj,quanj,βj,αj)
        h=min(hi,hj)
        if h==Inf
          qi=getQfromAsymptote(xi,βi,ci,αi,bi)
          qj=getQfromAsymptote(xj,βj,cj,αj,bj)
        else
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0.0
            println("Δ=0")
            h=h-1e-9
          end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ

        end
      else

        for i =1:3# 3 ord1 ,7 ord2
          acceptedi[i][1]=0.0; acceptedi[i][2]=0.0
          acceptedj[i][1]=0.0; acceptedj[i][2]=0.0
        end 
        for i=1:4
          cacheRootsi[i]=0.0
          cacheRootsj[i]=0.0
        end


      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))

   


      resi1,resi2=quadRootv2(coefi)
      resi3,resi4=quadRootv2(coefi2)
      resj1,resj2=quadRootv2(coefj)
      resj3,resj4=quadRootv2(coefj2)
   
      printTime=0.0
      if  resi1>0 && resi2>0 && resi3>0 && resi4>0
        @show αi,βi,quani,ci,bi,resi1,resi2,resi3,resi4
        printTime=simt
      end
      if  resj1>0 && resj2>0 && resj3>0 && resj4>0
        @show αj,βj,quanj,cj,bj,resj1,resj2,resj3,resj4
        printTime=simt
      end
  
      #construct intervals
      constructIntrval(acceptedi,resi1,resi2,resi3,resi4)

      constructIntrval(acceptedj,resj1,resj2,resj3,resj4)

      #find best H (largest overlap)
      ki=3;kj=3
      minF=0.0;segCounter=0
  
      h=0.0 #for testing... not needed
      while true
            currentLi=acceptedi[ki][1];currentHi=acceptedi[ki][2]
            currentLj=acceptedj[kj][1];currentHj=acceptedj[kj][2]
            if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                  b0=min(currentHi,currentHj )
                  a0=max(currentLj,currentLi)
                  if simt==printTime
                    @show a0,b0
                  end
                  if a0==0.0 && b0==Inf
                    qi=getQfromAsymptote(xi,βi,ci,αi,bi)
                    qj=getQfromAsymptote(xj,βj,cj,αj,bj)
                    #println("asymptote")
                    h=Inf
                    break # break because both acceptedintervals contain 1 interval
                  elseif a0==0.0 && b0!=Inf
                    
                    h0,tempphi=goldenSearch(a0,b0,50,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                    temp2[1]+=1
                    if tempphi<minF 
                      if minF!=0.0 # this is not needed since only here at first iteration
                        temp3[1]+=1
                      end
                      minF=tempphi
                      h=h0
                      if simt==printTime
                        @show h
                      end
                      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                      if Δ==0.0
                        println("Δ=0")
                        h=h-1e-9
                      end
                        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                    end
                  
                   
                    break # break because both intervals reached 0 at left
                  elseif a0!=0.0 && b0==Inf
                    h0,tempphi=goldenSearch(a0,a0*50.0,50,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                    temp2[1]+=1
                    if tempphi<minF  # this is not needed since only here at first iteration
                      if minF!=0.0 # this is not needed since only here at first iteration
                        temp3[1]+=1
                        end
                      minF=tempphi
                      h=h0
                      if simt==printTime
                        @show h
                      end
                      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                      if Δ==0.0
                        println("Δ=0")
                        h=h-1e-9
                      end
                        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                    end
                        if currentLj<currentLi#remove last seg of resi
                          ki-=1
                        else #remove last seg of resj
                          kj-=1
                        end
                       # println("1st iter a0,inf simt=$simt minf=$minF cntr=$segCounter")
                        segCounter+=1
                  elseif a0!=0.0 && b0!=Inf
                   #=  if simt==printTime
                      @show a0,b0
                    end =#
                    h0,tempphi=goldenSearch(a0,b0,50,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                    temp2[1]+=1
                    if tempphi<minF  # this is not needed since only here at first iteration
                      if minF!=0.0
                        temp3[1]+=1
                      end
                      minF=tempphi
                      h=h0
                      if simt==printTime
                        @show h
                      end
                      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                      if Δ==0.0
                        println("Δ=0")
                        h=h-1e-9
                      end
                        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                    end
                    if currentLj<currentLi#remove last seg of resi
                      ki-=1
                    else #remove last seg of resj
                      kj-=1
                    end
                  end
                 # println("comarebounds simt=$simt minf=$minF cntr=$segCounter"); 
                  segCounter+=1          
            else
                      if currentHj==0.0 && currentHi==0.0#empty
                        ki-=1;kj-=1
                      elseif currentHj==0.0 && currentHi!=0.0
                        kj-=1
                      elseif currentHi==0.0 && currentHj!=0.0
                        ki-=1
                      else #both non zero
                            if currentLj<currentLi#remove last seg of resi
                              ki-=1
                            else #remove last seg of resj
                              kj-=1
                            end
                      end
            end

      end # end while
    end#end ifelse bc>0
        #=   if simt==printTime
            @show "fin",h

          end =#
      q[j][0]=qj
      tq[j]=simt 
      q[index][0]=qi# store back helper vars
      trackSimul[1]+=1 # do not have to recomputeNext if qi never changed
     
  
  end
  return iscycle
end   
 #phi(xi::Float64,qi::Float64,xj::Float64,qj::Float64)=-((xi-qi)^2+(xj-qj)^2)

 function phi(h::Float64,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
  if Δ==0.0
    h=h-1e-9
  end
    qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
    qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
    return -((abs(xi-qi)/quani)+(abs(xj-qj)/quanj))
 
 #=  if Δ!=0.0
    qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
    qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
    return -((abs(xi-qi)/quani)+(abs(xj-qj)/quanj))
  else
    println("phi returns 0")
    return 0.0
  end =#

 end
 #= function goldenSearch(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  ρ=0.382
  a1=a0;b1=b0
  for k=1:N
    a0=a1;b0=b1
    a1=a1+ρ*(b1-a1);b1=b1-ρ*(b1-a1);
    if phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)<phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
      a1=a0;
    else
      b1=b0
    end
  end
  return (a1+b1)/2
 end =#


 #special golden search...if it does not see an improvment it stops with  previous option 
#=  function goldenSearch0(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  ρ=0.382
  a1=a0;b1=b0
 
  counterLeft=0;counterRight=0;
  for k=1:N
    a0=a1;b0=b1
    a1=a1+ρ*(b1-a1);b1=b1-ρ*(b1-a1);
    phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
    if phia<phib
      a1=a0;counterLeft+=1
    else
      b1=b0;counterRight+=1
    end
    #= if counterLeft>counterRight+20 # 
      println("return a0")
      return a0,phia
    elseif counterLeft+20<counterRight
      return b0,phib
    end =#
  end
  phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
  if phia<phib
   # @show a0,a1,b0,b1
    return a1,phia
  else
    return b1,phib
  end
  #println("return (a1+b1)/2")
 
 end =#
 function goldenSearch(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  ρ=0.382
  a1=a0;b1=b0
  counterLeft=0;counterRight=0;
  for k=1:N
    a0=a1;b0=b1
    a1=a1+ρ*(b1-a1);b1=b1-ρ*(b1-a1);
    phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
    if phia<phib
      a1=a0;counterLeft+=1
    else
      b1=b0;counterRight+=1
    end
    #= if counterLeft>counterRight+20 # 
      println("return a0")
      return a0,phia
    elseif counterLeft+20<counterRight
      return b0,phib
    end =#
  end
  phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
  if phia<phib
   # @show a1,a0
    return a1,phia
  else
    return b1,phib
  end
  #println("return (a1+b1)/2")
 
 end
#=  function goldenSearch1(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  ρ=0.382
  a1=a0;b1=b0
  counterLeft=0;counterRight=0;
  for k=1:N
    a0=a1;b0=b1
    a1=a1+ρ*(b1-a1);b1=b1-ρ*(b1-a1);
    phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
    if phia<phib
      a1=a0;counterLeft+=1
    else
      b1=b0;counterRight+=1
    end
    #= if counterLeft>counterRight+20 # 
      println("return a0")
      return a0,phia
    elseif counterLeft+20<counterRight
      return b0,phib
    end =#
  end
  phia=phi(a1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b1,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
  if phia<phib
    
    return a1,phia
  else
   # @show b1
    return b1,phib
  end
  #println("return (a1+b1)/2")
 
 end =#
 function goldenSearch2(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  phia=phi(a0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
  #@show a0,b0,phia,phib
  if phia<phib
      println("return *******************************************a0")
      return a0,phia

    else
      return b0,phib
    end
 end 

function nmisCycle_and_simulUpdate(cacheRootsi::Vector{Float64},cacheRootsj::Vector{Float64},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,rejected,simpleCase,temp2,temp3)
  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index,simt)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,d,cacheA,j,j,simt)
  ajj=cacheA[1]
  α=-aii-ajj
    β=aii*ajj-aij*aji
    if α<0 || β<0
      rejected[1]+=1
      return false
    end

   xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  #qaux[j][1]=qj;
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;

  xj=x[j][0]
  tx[j]=simt

 qiminus=qaux[index][1]
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qiminus
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qiminus-aij*qj
    iscycle=false
    dxj=aji*qi+ajj*qj+uji #only future qi   #emulate fj
   
    dxithrow=aii*qi+aij*qj+uij #only future qi
                                                      
  qjplus=xj+sign(dxj)*quanj  #emulate updateQ(j)...

    dxi=aii*qi+aij*qjplus+uij #both future qi & qj   #emulate fi
  
  ########old condition:Union 
     #= if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
    if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
        iscycle=true
    end
  end   =#

    ########condition:Union i union
    if (abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0)
      cancelCriteria=1e-6*quani
      if abs(dxi)>3*abs(dxithrow) || abs(dxi)*3<abs(dxithrow) ||  (dxi*dxithrow)<0.0  
          iscycle=true
          if abs(dxi)<cancelCriteria && abs(dxithrow)<cancelCriteria #significant change is relative to der values. cancel if all values are small
            iscycle=false
           # println("cancel dxi & dxthr small")
          end
      end
      if  abs(dxi)>10*abs(ẋi) || abs(dxi)*10<abs(ẋi) #||  (dxi*ẋi)<0.0 
        iscycle=true
        if abs(dxi)<cancelCriteria && abs(ẋi)<cancelCriteria #significant change is relative to der values. cancel if all values are small
          iscycle=false
         # println("cancel dxi & ẋi small")
        end
      end
      
      if abs(dxj)<1e-6*quanj && abs(ẋj)<1e-6*quanj #significant change is relative to der values. cancel if all values are small
        iscycle=false
        #println("cancel dxj & ẋj small")
      end

    end   
   

   # iscycle=false
  if iscycle

   
      
       #find positive zeros f=+-Δ
      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij
    

      if bi*ci>0.0 && bj*cj>0.0
        simpleCase[1]+=1
        hi=bcPOS(bi,ci,quani,βi,αi)
        hj=bcPOS(bj,cj,quanj,βj,αj)
        h=min(hi,hj)
        if h==Inf
          qi=getQfromAsymptote(xi,βi,ci,αi,bi)
          qj=getQfromAsymptote(xj,βj,cj,αj,bj)
        else
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0.0
            println("Δ=0")
            h=h-1e-9
          end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ

        end
      else

        for i =1:3# 3 ord1 ,7 ord2
          acceptedi[i][1]=0.0; acceptedi[i][2]=0.0
          acceptedj[i][1]=0.0; acceptedj[i][2]=0.0
        end 
        for i=1:4
          cacheRootsi[i]=0.0
          cacheRootsj[i]=0.0
        end
       


      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))

   


      resi1,resi2=quadRootv2(coefi)
      resi3,resi4=quadRootv2(coefi2)
      resj1,resj2=quadRootv2(coefj)
      resj3,resj4=quadRootv2(coefj2)
   
      printTime=0.0
      if  resi1>0 && resi2>0 && resi3>0 && resi4>0
        @show αi,βi,quani,ci,bi,resi1,resi2,resi3,resi4
        printTime=simt
      end
      if  resj1>0 && resj2>0 && resj3>0 && resj4>0
        @show αj,βj,quanj,cj,bj,resj1,resj2,resj3,resj4
        printTime=simt
      end
  
      #construct intervals
      constructIntrval(acceptedi,resi1,resi2,resi3,resi4)

      constructIntrval(acceptedj,resj1,resj2,resj3,resj4)

      #find best H (largest overlap)
      ki=3;kj=3
      minF=0.0;segCounter=0
  
      h=0.0 #for testing... not needed
      while true
            currentLi=acceptedi[ki][1];currentHi=acceptedi[ki][2]
            currentLj=acceptedj[kj][1];currentHj=acceptedj[kj][2]
            if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                  b0=min(currentHi,currentHj )
                  a0=max(currentLj,currentLi)
                  if simt==printTime
                    @show a0,b0
                  end
                  if a0==0.0 && b0==Inf
                    qi=getQfromAsymptote(xi,βi,ci,αi,bi)
                    qj=getQfromAsymptote(xj,βj,cj,αj,bj)
                    #println("asymptote")
                    h=Inf
                    break # break because both acceptedintervals contain 1 interval
                  elseif a0==0.0 && b0!=Inf
                    tempphi=phi(b0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                   
                    if tempphi<minF
                     
                      h=b0; 
                      if simt==printTime#6.020028519408029 
                        @show h
                      end
                      if minF!=0
                        println("minF!=0")
                      end
                      if segCounter!=0  
                       # println("reached 0 at left simt=$simt tempphi=$tempphi cntr=$segCounter") 
                      end
                          minF=tempphi
                          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                          if Δ==0.0
                            println("Δ=0")
                            h=h-1e-9
                          end
                            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      
                          
                    end
                  
                   
                    break # break because both intervals reached 0 at left
                  elseif a0!=0.0 && b0==Inf
                        tempphi=phi(a0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                        
                        if tempphi<minF  # this is not needed since only here at first iteration
                         
                          h=a0;
                          minF=tempphi
                          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                          if Δ==0.0
                            println("Δ=0")
                            h=h-1e-9
                          end
                            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                        end
                        if currentLj<currentLi#remove last seg of resi
                          ki-=1
                        else #remove last seg of resj
                          kj-=1
                        end
                       # println("1st iter a0,inf simt=$simt minf=$minF cntr=$segCounter")
                        segCounter+=1
                  elseif a0!=0.0 && b0!=Inf
                    if simt==printTime
                      @show a0,b0
                    end
                    h0,tempphi=compareBounds(a0,b0,50,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
                    temp2[1]+=1
                    if tempphi<minF  # this is not needed since only here at first iteration
                      if minF!=0.0 # this is not needed since only here at first iteration
                        temp3[1]+=1
                        end
                      minF=tempphi
                      h=h0
                      if simt==printTime
                        @show h
                      end
                      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                      if Δ==0.0
                        println("Δ=0")
                        h=h-1e-9
                      end
                        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                    end
                    if currentLj<currentLi#remove last seg of resi
                      ki-=1
                    else #remove last seg of resj
                      kj-=1
                    end
                  end
                 # println("comarebounds simt=$simt minf=$minF cntr=$segCounter"); 
                  segCounter+=1          
            else
                      if currentHj==0.0 && currentHi==0.0#empty
                        ki-=1;kj-=1
                      elseif currentHj==0.0 && currentHi!=0.0
                        kj-=1
                      elseif currentHi==0.0 && currentHj!=0.0
                        ki-=1
                      else #both non zero
                            if currentLj<currentLi#remove last seg of resi
                              ki-=1
                            else #remove last seg of resj
                              kj-=1
                            end
                      end
            end

      end # end while
    end#end ifelse bc>0
        #=   if simt==printTime
            @show "fin",h

          end =#
      q[j][0]=qj
      tq[j]=simt 
      q[index][0]=qi# store back helper vars
      trackSimul[1]+=1 # do not have to recomputeNext if qi never changed
     
  
  end
  return iscycle
end   
 #not search just compare bounds
 function compareBounds(a0::Float64,b0::Float64,N::Int,aii::Float64,ajj::Float64,aij::Float64,aji::Float64,uij::Float64,uji::Float64,xi::Float64,xj::Float64,quani::Float64,quanj::Float64)
  phia=phi(a0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj);phib=phi(b0,aii,ajj,aij,aji,uij,uji,xi,xj,quani,quanj)
  #@show a0,b0,phia,phib
  if phia<phib
      println("return *******************************************a0")
      return a0,phia

    else
      return b0,phib
    end
 end 

function bcPOS(b,c,quan,β,α)
  if c>0
    if quan<c/β
      coef=NTuple{3,Float64}((β*quan-c,α*quan-b,quan))
      h=minPosRootv1(coef)
    else
      h=Inf
    end
  else
    if -quan>c/β
      coef=NTuple{3,Float64}((-β*quan-c,-α*quan-b,-quan))
      h=minPosRootv1(coef)
    else
      h=Inf
    end
  end
return h

end
 function hasThreePos(a,b,c,d)
  if a>0 && b>0 && c>0 
    return true
  end
  if a>0 && b>0 && d>0 
    return true
  end
  if a>0 && c>0 && d>0 
    return true
  end
  if b>0 && c>0 && d>0 
    return true
  end
  return false
 end
