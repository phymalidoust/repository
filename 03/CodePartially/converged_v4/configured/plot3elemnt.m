

load cpu_0.txt; load cpu_1.txt; load cpu_2.txt; load cpu_3.txt
load cpu_4.txt; load cpu_5.txt; load cpu_6.txt; load cpu_7.txt
load cpu_8.txt; load cpu_9.txt; load cpu_10.txt; load cpu_11.txt
load cpu_12.txt; load cpu_13.txt; load cpu_14.txt



MCEnrgy = [cpu_0;cpu_1;cpu_2;cpu_3;
          cpu_4;cpu_5;cpu_6;cpu_7;
          cpu_8;cpu_9;cpu_10;cpu_11;
          cpu_12;cpu_13;cpu_14];

% alvc = [0.001:0.0225:1.0];
% alvc(1) = 0.001;
% alvc(length(alvc)) = 0.999;

alvc = [0.001:0.0225:1.0];
alvc(1) = 0.003;
alvc(length(alvc)) = 0.997;

nn = 0;

for nal = 1:length(alvc)
    
    al0 = alvc(nal);
    
    znvc = [0.001:0.0225:1.0 - al0 ];
    znvc2 = [0.001:0.001:1.0 - al0 ];
    znvc(length(znvc)) = znvc2(length(znvc2));
    znvc(1) = 0.003;
    
    for nzn = 1:length(znvc)
        zn0 = znvc(nzn);
        mg0 = 1.0 - al0 - zn0;
        %mgvc(nzn) = mg0;
        nn = nn + 1;
%    res(nal, nzn) = MCEnrgy(nn,1,1,1) - al0*EAl - zn0*EZn - mg0*EMg;
ordEngy(nal, nzn) = MCEnrgy(nn,1);
% res(nal, nzn) = MCEnrgy(nn,1) - MCEnrgy(nn,2)*EAl - MCEnrgy(nn,3)*EZn - MCEnrgy(nn,4)*EMg;

alum(nal, nzn) =  MCEnrgy(nn,2);
zink(nal, nzn) =  MCEnrgy(nn,3);
magn(nal, nzn) =  MCEnrgy(nn,4);

    end
    
end

% X = alum;
% Y = zink;
% Z = magn;

% surf(ordEngy)
% surf(X,Y,ordEngy)
% surf(Y,X,ordEngy)
% surf(Z,Y,ordEngy)
% surf(X,Z,ordEngy)
% surf(Z,X,ordEngy)

% alvc = [0.05:0.0207:0.92];
% for jj=1:length(alvc)
% X(:,jj) = alvc';
% end
% for jj=1:length(alvc)
% Y(jj,:) = alvc;
% end
% 
% alvc = [0, 0.05:0.0207:0.92, 1.0];
% for jj=1:length(alvc)
% Xq(:,jj) = alvc';
% end
% for jj=1:length(alvc)
% Yq(jj,:) = alvc;
% end
% 
% Vq = interp2(X,Y,ordEngy,Xq,Yq,'linear','extrap');
% 
% Vq=interpextrap2(X,Y,ordEngy,Xq,Yq);

% 
% Vq = interp2(ordEngy,Xq,Yq);

% x = [mang(1,2),mang(1,1)];
% 
% for kk=1:length(alvc)-1
% 
% v = [ordEngy(kk,2),ordEngy(kk,1)];
% 
% xq = [x,1];
% vq = interp1(x,v,xq,'linear','extrap');
% Elmn(kk) = vq(3);
% 
% end
% 
% xA = [0,Energy(:,2)',1];
% Energy_ex = [Elm1,Energy(:,1)',Elmn];
% 
% 
% al0 = alvc(1);
% znvc = [0.05:0.0207:1.0 - al0 - 0.01];
% EMg1 = ordEngy(1,1);
% EZn1 = ordEngy(1,length(znvc));
% EAl1 = ordEngy(length(alvc)-1,1);
% 
% 
% matrix = [0.90 , 0.05 , 0.05;
%           0.015 , 0.955 , 0.03;
%           0.015 , 0.03 , 0.9194];
% z1 = [EMg1 ; EZn1 ; EAl1];      
% z2 = inv(matrix)*z1;     
% 
% EMg = z2(1);
% EZn = z2(2);
% EAl = z2(3);




EMg = ordEngy(1,1);
EZn = ordEngy(1,length(ordEngy(1,:)));
EAl = ordEngy(length(ordEngy(:,1)),1);

nn = 0;

for nal = 1:length(alvc)
    
    al0 = alvc(nal);
    
    znvc = [0.001:0.0225:1.0 - al0 ];
    znvc2 = [0.001:0.001:1.0 - al0 ];
    znvc(length(znvc)) = znvc2(length(znvc2));
    znvc(1) = 0.003;
    
    for nzn = 1:length(znvc)
        zn0 = znvc(nzn);
        mg0 = 1 - al0 - zn0;
        %mgvc(nzn) = mg0;
        nn = nn + 1;
%    res(nal, nzn) = MCEnrgy(nn,1,1,1) - al0*EAl - zn0*EZn - mg0*EMg;
%ordEngy(nal, nzn) = MCEnrgy(nn,1);
res(nal, nzn) = MCEnrgy(nn,1) - MCEnrgy(nn,2)*EAl - MCEnrgy(nn,3)*EZn - MCEnrgy(nn,4)*EMg;

alum(nal, nzn) =  MCEnrgy(nn,2);
zink(nal, nzn) =  MCEnrgy(nn,3);
magn(nal, nzn) =  MCEnrgy(nn,4);

    end
    
end

figure
surf(res./1000)
% figure
% surf(ordEngy./2000)
figure
contourf(res./1000,15)

% clear all

figure
for aa = 1:length(alvc)
    hold all; 
    subplot(221);plot(nonzeros(alum(:,aa)),nonzeros(res(:,aa))./2000);xlabel('X_{Al}%')
    hold all; 
    subplot(223);plot(nonzeros(magn(:,aa)),nonzeros(res(:,aa))./2000);xlabel('X_{Mg}%')
    hold all; 
    subplot(222);plot(nonzeros(zink(aa,:)),nonzeros(res(aa,:))./2000);xlabel('X_{Zn}%')
    hold all; 
    subplot(224);plot(nonzeros(magn(aa,:)),nonzeros(res(aa,:))./2000);xlabel('X_{Mg}%')

    pause
end





% alvc = [0.01:0.0165:1.0];
% nn = 0;
% 
% for nal = 1:length(alvc)
%     
%     al0 = alvc(nal);
%     znvc = [0.01:0.0165:1.0 - al0 ];
%     
%     for nzn = 1:length(znvc)
%         zn0 = znvc(nzn);
%         mg0 = 1.0 - al0 - zn0;
%         %mgvc(nzn) = mg0;
%         nn = nn + 1;
% 
% alum(nal, nzn) =  al0;
% zink(nal, nzn) =  zn0;
% magn(nal, nzn) =  zn0;
% 
%     end
%     
% end

clear all
