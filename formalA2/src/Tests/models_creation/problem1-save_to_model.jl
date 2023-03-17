using formalA

function test()
    odeprob = @NLodeProblem begin
        
      u = [1.0, 0.0,0.0]
      discrete = [0.0]
      du[1] = 100.8*(   9.523809523809524e-5 * u[2] - u[1] * u[2] + u[1] * ( 1.0 - u[1] ) ) 
      du[2] = 40320.0*( - 9.523809523809524e-5 * u[2] - u[1] * u[2] +  u[3] ) 
      du[3] = u[1] - u[3]
    end
  display(odeprob)
    save_prob_to_model(odeprob,"/home/unknown/relaxedqssA/relaxedqssA/src/models/NonLinear.jl","oregonator") #any location you want
   
end
test()
