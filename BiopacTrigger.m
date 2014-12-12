function trigger = BiopacTrigger(emg, slope)

demg = diff(emg);
demg = abs(demg);
[x y] = find(demg > slope)
dx = diff(x)

if length(x) == 0
    trigger = [];
else

C{1} = x(1);
j = 1;
cont = 1;
mx = [];

for i = 1:length(dx)
    
    if dx(i) < 100
        C{j} = [C{j}, x(i+1)];
    else
        j = j+1;
        C{j} = x(i+1);
    end
end

clear i

for i = 1:length(C)
    m(i) = max(C{i});
end

clear i


for i = 1:length(m)
    absy(i) = max(abs(emg(m(i)-20:m(i)+20)));
    [X Y] = find(abs(emg) == absy(i));
    
    if length(X)>1
        
        dX = diff(X);
        cont = 2;
        aux = X(1);
        for w = 1:length(dX)
            if dX(w)>100
                aux(cont) = X(w+1);
                cont = cont+1;
            end
                
        end
        
        X = aux';
      
        
        for j = 1:length(X)
            
            if ismember(X(j),mx) == 0
                mx(i) = X(j);
                break
            end
        end
        
    else
        mx(i) = X;
    end
    
end



for i = 1:length(mx)
    my(i) = emg(mx(i));
end

trigger = [mx' my'];
end

% figure
% plot(emg,'.-')
% hold on
% 
% plot(trigger(:,1),trigger(:,2),'ro')