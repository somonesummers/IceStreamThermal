function tau_base = channelSpacing(y,lambda,w,amp)
%Create strength distribution for channels spaced (lambda) apart, with width
%(w), and amplitude (amp)

channel =@(loc,amp,w) amp*(heaviside(y-loc+0.5*w).*...
                           cos(pi*(y-loc)/w).^2.*...
                           heaviside(loc+0.5*w-y));
                       
tau_base = zeros(length(y),1);
for i = lambda:lambda:max(y)
    tau_base = tau_base + channel(i,amp,w);
end

end

