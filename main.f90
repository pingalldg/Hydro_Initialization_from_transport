real*8::tau_ampt,dx,detas,etas_min,etas_wind,sigma,sigma_etas,max_etas,ks_0
integer::k,netas,netas_print,input_method
open(777,file="input_file.dat")
read(777,*)input_method,k,dx,tau_ampt,sigma,netas,detas,etas_min,etas_wind,sigma_etas,max_etas,ks_0
if (input_method.eq.1)then
call ampttohydro1(k,dx,tau_ampt,sigma,1,netas,detas,etas_min,etas_wind,ks_0)
elseif(input_method.eq.2)then
call ampttohydro1(k,dx,tau_ampt,sigma,netas,netas,detas,etas_min,etas_wind,ks_0)
elseif(input_method.eq.3)then
call ampttohydro2(k,dx,tau_ampt,sigma,netas,netas,detas,etas_min,-max_etas,max_etas,sigma_etas,ks_0)
end if
 close(777)
end
