        function imf = dsp_fin(x)

c = x(:)'; % copy of the input signal (as a row vector)
N = length(x);

%-------------------------------------------------------------------------
% loop to decompose the input signal into successive IMF

imf = []; % Matrix which will contain the successive IMF, and the residue

while (1) % the stop criterion is tested at the end of the loop
   
   %-------------------------------------------------------------------------
   % inner loop to find each imf
   
   h = c; % at the beginning of the sifting process, h is the signal
   SD = 1; % Standard deviation which will be used to stop the sifting process
   
   while SD > 0.3
      % while the standard deviation is higher than 0.3 (typical value)
      
      % find local max/min points
      d = diff(h); % approximate derivative
      maxmin = []; % to store the optima (min and max without distinction so far)
      for i=1:N-2
         if d(i)==0                        % we are on a zero
            maxmin = [maxmin, i];
         elseif sign(d(i))~=sign(d(i+1))   % we are straddling a zero so
            maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
         end
      end
      
      if size(maxmin,2) < 2 % then it is the residue
         break
      end
      
      % divide maxmin into maxes and mins
      if maxmin(1)>maxmin(2)              % first one is a max not a min
         maxes = maxmin(1:2:length(maxmin));
         mins  = maxmin(2:2:length(maxmin));
      else                                % is the other way around
         maxes = maxmin(2:2:length(maxmin));
         mins  = maxmin(1:2:length(maxmin));
      end
      
      % make endpoints both maxes and mins
      maxes = [1 maxes N];
      mins  = [1 mins  N];
      
      
      %-------------------------------------------------------------------------
      % spline interpolate to get max and min envelopes; form imf
      maxenv = spline(maxes,h(maxes),1:N);
      minenv = spline(mins, h(mins),1:N);
      
      m = (maxenv + minenv)/2; % mean of max and min enveloppes
      prevh = h; % copy of the previous value of h before modifying it
      h = h - m; % substract mean to h
      
      % calculate standard deviation
      eps = 0.0000001; % to avoid zero values
      SD = sum ( ((prevh - h).^2) ./ (prevh.^2 + eps) );
      
   end
   
   imf = [imf; h]; % store the extracted IMF in the matrix imf
   % if size(maxmin,2)<2, then h is the residue
   %display(imf);
   % stop criterion of the algo.
   if size(maxmin,2) < 2
      break
   end
   
   c = c - h; % substract the extracted IMF from the signal
   
end

if( rem(size(imf,1), 2) == 0)
    numPlot = size(imf,1)/2;
else
    numPlot = (size(imf,1)+1)/2;
end

for i = 1:size(imf,1)
    subplot(numPlot,2,i);
    plot(imf(i,:));
end
fs = 1000;
t = 0:2/fs:2-1/fs;

energy = zeros(1,size(imf,1));
for i= 1:size(imf,1)
     %syms x y z energy pow;
     f=imf(i,:);
     energy(i) = f*f';
end
B = sort(energy);
tot_eng=0.4*(sum(energy));
temp = 0;
for i = 1:size(energy,2)
    temp = temp + B(i);
    if(temp<tot_eng)
        for j=1:size(energy,2)
            if(energy(j) == B(i))
                z = hilbert(imf(j,:));
                instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
                figure;
                plot(smooth(instfreq));
                xlabel('Time')
                ylabel('Hz')
                ylim([0 500]);
                grid on
                title('Instantaneous Frequency')
            end
        end
    else
        break
    end
end
org=sum(imf);

figure;
subplot(1,2,1);
plot(org);
title('sum of IMFs');
subplot(1,2,2);
plot(x);
title('original');
return