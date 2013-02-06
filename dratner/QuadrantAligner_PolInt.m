function QuadrantAligner_PolInt(img)
%-------------------------------------------------------------
% QuadrantAligner
% Finds alignment of quadrants in CSPAD from a powder pattern
%
% D. Ratner, May 13, 2012
%-------------------------------------------------------------

myq=440:.5:490;     % q range of peak in pixels (Cu run14)
%myq=250:.5:320;     % q range of peak in pixels (Ti run96, ring1)
myq=330:.25:365;     % q range of peak in pixels (Ti run96, ring2)
%myq = 535:.5:575;



Nint=2;             % Interpolation factor (increase points by Nint^2)

qD=1;               % Range of center to search (look for center +/- qD)
qd=.25;             % Delta step for center search (-qD:qd:qD)


% read in binary file.
%filename='xtcexplorer-r0097_image_ev1.dat'; binFile = fopen(filename,'r');
%temp_img2 = fread(binFile,[1800,1800],'short');


% % estimated centers of each quadrant (Cu run14)
% center1 = [850.5, 854];
% center2 = [855.5, 852.5];
% center3 = [852, 858];
% center4 = [858.5, 857.5];

% % estimated centers of each quadrant (Cu run14)
% center1 = [875.5, 872];
% center2 = [883, 872.5];
% center3 = [880, 876];
% center4 = [883, 878];


% estimated centers of each quadrant (Ti run96, ring 2)
center1 = [877, 872.5];
center2 = [882.5, 873.5];
center3 = [879, 878];
center4 = [883, 877.5];

%center=[881 876]; center1 = center; center2 = center; center3 = center; center4 = center;


% % make test ring;
% myq=140:.3:190;
% center=[857 855]; r=170; delta_r=.5;
% img=1e3*test_aliasing(img,center,r,delta_r); 
% center1 = center; center2 = center; center3 = center; center4 = center;

[nx,ny]=size(img);
[X,Y] = meshgrid(1:nx,1:ny);

% dividing line of quadrants
q1=912; q2=850; 
%q1=1799; q2=1799;

% separate out quadrants
quad1 = img(1:q1,1:q2);  X1=X(1:q1,1:q2); Y1=Y(1:q1,1:q2);
quad2 = img(1:q2,q2+1:end); X2=X(1:q2,q2+1:end); Y2=Y(1:q2,q2+1:end);
quad3 = img(q1+1:end,1:q1); X3=X(q1+1:end,1:q1); Y3=Y(q1+1:end,1:q1);
quad4 = img(q2+1:end,q1+1:end); X4=X(q2+1:end,q1+1:end); Y4=Y(q2+1:end,q1+1:end);

% replot quadrants according to center1-4
replot_quad(quad1,quad2,quad3,quad4,center1,center2,center3,center4,X1,X2,X3,X4,Y1,Y2,Y3,Y4);

disp('Starting search for quadrant centers');
drawnow

% search for center of each quadratn
tic; [A1,W1,Q1]=find_center(quad1,X1,Y1,center1,Nint,myq,qd,qD); disp(['quad1: elapsed time ' num2str(toc) ' seconds']);
tic; [A2,W2,Q2]=find_center(quad2,X2,Y2,center2,Nint,myq,qd,qD); disp(['quad2: elapsed time ' num2str(toc) ' seconds']);
tic; [A3,W3,Q3]=find_center(quad3,X3,Y3,center3,Nint,myq,qd,qD); disp(['quad3: elapsed time ' num2str(toc) ' seconds']);
tic; [A4,W4,Q4]=find_center(quad4,X4,Y4,center4,Nint,myq,qd,qD); disp(['quad4: elapsed time ' num2str(toc) ' seconds']);

% vectors for center search
XC1=center1(1)-qD:qd:center1(1)+qD; YC1=center1(2)-qD:qd:center1(2)+qD;
XC2=center2(1)-qD:qd:center2(1)+qD; YC2=center2(2)-qD:qd:center2(2)+qD;
XC3=center3(1)-qD:qd:center3(1)+qD; YC3=center3(2)-qD:qd:center3(2)+qD;
XC4=center4(1)-qD:qd:center4(1)+qD; YC4=center4(2)-qD:qd:center4(2)+qD;

% figure(1); imagesc(XC1,YC1,A1); title('Quadrant 1');
% figure(2); imagesc(XC2,YC2,A2); title('Quadrant 2');
% figure(3); imagesc(XC3,YC3,A3); title('Quadrant 3');
% figure(4); imagesc(XC4,YC4,A4); title('Quadrant 4');
% 
% figure(11); imagesc(XC1,YC1,W1); title('Quadrant 1');
% figure(12); imagesc(XC2,YC2,W2); title('Quadrant 2');
% figure(13); imagesc(XC3,YC3,W3); title('Quadrant 3');
% figure(14); imagesc(XC4,YC4,W4); title('Quadrant 4');

% Plot peak value of ring
figure(21); 
subplot(2,2,1); imagesc(XC1,YC1,A1); title('Quadrant 1, Ring peak'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,2); imagesc(XC2,YC2,A2); title('Quadrant 2, Ring peak'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,3); imagesc(XC3,YC3,A3); title('Quadrant 3, Ring peak'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,4); imagesc(XC4,YC4,A4); title('Quadrant 4, Ring peak'); xlabel('X (pixels)'); ylabel('Y (pixels)')

% Plot width of ring
figure(22);
subplot(2,2,1); imagesc(XC1,YC1,W1); title('Quadrant 1, Ring width');
subplot(2,2,2); imagesc(XC2,YC2,W2); title('Quadrant 2, Ring width');
subplot(2,2,3); imagesc(XC3,YC3,W3); title('Quadrant 3, Ring width');
subplot(2,2,4); imagesc(XC4,YC4,W4); title('Quadrant 4, Ring width');

% Plot radius of ring
figure(23); 
subplot(2,2,1); imagesc(XC1,YC1,Q1); title('Quadrant 1, Ring radius'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,2); imagesc(XC2,YC2,Q2); title('Quadrant 2, Ring radius'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,3); imagesc(XC3,YC3,Q3); title('Quadrant 3, Ring radius'); xlabel('X (pixels)'); ylabel('Y (pixels)')
subplot(2,2,4); imagesc(XC4,YC4,Q4); title('Quadrant 4, Ring radius'); xlabel('X (pixels)'); ylabel('Y (pixels)')



function [A,W,Q]=find_center(img,X,Y,center,Nint,myq,qd,qD)

% vectors of centers to try
XC=center(1)-qD:qd:center(1)+qD;
YC=center(2)-qD:qd:center(2)+qD;

% pick out img pixels within range myq
Rtemp=sqrt((X-center(1)).^2+(Y-center(2)).^2);      % convert to polar
Itemp(1,:)=Rtemp(:)'; Itemp(2,:)=X(:)'; Itemp(3,:)=Y(:)'; Itemp(4,:)=img(:)';  % make vector structure 
Itemp = sortrows(Itemp')';      % sort by radius
kmin = find(Itemp(1,:) > myq(1),1);     % find ring minimum
kmax = find(Itemp(1,:) > myq(end),1);   % find ring maximum
X_ring=Itemp(2,kmin:kmax); Y_ring=Itemp(3,kmin:kmax);  img_ring=Itemp(4,kmin:kmax);

% make finer grid for interpolation
[nx,ny]=size(img);
x_f=linspace(min(min(X)),max(max(X)),nx*Nint);
y_f=linspace(min(min(Y)),max(max(Y)),ny*Nint);
[Xf,Yf]=meshgrid(x_f,y_f);

% select pixels from fine grid within range myq
Rtemp=sqrt((Xf-center(1)).^2+(Yf-center(2)).^2);
clear('Itemp');
Itemp(1,:)=Rtemp(:)'; Itemp(2,:)=Xf(:)'; Itemp(3,:)=Yf(:)';  % make vector structure 
Itemp = sortrows(Itemp')';      % sort by radius
kmin = find(Itemp(1,:) > myq(1),1);     % find ring minimum
kmax = find(Itemp(1,:) > myq(end),1);   % find ring maximum
XI=Itemp(2,kmin:kmax); YI=Itemp(3,kmin:kmax); 

% interpolate img within finer grid
F = TriScatteredInterp(X_ring',Y_ring',img_ring','natural');  % interpolation structure
img_f = F(XI',YI');        % evaluate structure on fine grid

max_A=0; min_W=1e3;
W=zeros(length(YC),length(XC)); A=zeros(length(YC),length(XC)); Q=A;
for l=1:length(XC)
    for p=1:length(YC)
        
        % find radius of each pixel
        R=sqrt((XI-XC(l)).^2+(YI-YC(p)).^2);

        % convert to vector and sort by radius
        Ivec(1,:)=R; Ivec(2,:)=img_f; 
        Ivec=sortrows(Ivec')';

        % find average pixel value for each q bin (myq)
        max_k = length(Ivec);        
        k=find(Ivec(1,:) > myq(1),1);   
        myproj = zeros(size(myq));
        for j=1:length(myq)
            n=0;
            while k<max_k && Ivec(1,k) < myq(j)
                myproj(j) = myproj(j)+Ivec(2,k);
                n=n+1;
                k=k+1;
            end
            if n>0
                myproj(j)=myproj(j)/n;  % normalize by number of pixels in ring
            end
        end
        
        % find maximum q value
        [A(p,l),imax]=max(myproj);        

        %figure(109); hold on; plot(myq,myproj)
        
        % find width of peak
        edge1=imax;
        while edge1>0 && myproj(edge1)>A(p,l)/2
            edge1=edge1-1;
        end        
        edge2=imax;
        while edge2<=length(myproj) && myproj(edge2)>A(p,l)/2
            edge2=edge2+1;
        end
        W(p,l)=edge2-edge1;

        % record center pixel with maximum peak amplitude
        if A(p,l)>max_A            
            maxmax=imax;
            max_A=A(p,l);
        end

        % record center pixel with minimum peak width
        if W(p,l)<min_W            
            minmin=imax;
            min_W=W(p,l);
        end        
        
        Q(p,l) = myq(imax);
        
    end
end
disp(['peak value of ' num2str(max_A) ' at Q=' num2str(myq(maxmax)) ', minimum width of ' num2str(min_W) ' at Q=' num2str(myq(minmin))]);


function test_img=test_aliasing(img,center,r,r_sig)
% make gaussian ring for testing code

xc=center(1,1);
yc=center(1,2);

[m,n]=size(img);
x=(1:m)-xc;
y=(1:n)-yc;

[X,Y]=meshgrid(x,y);

R=sqrt(X.^2+Y.^2);

test_img=exp(-(R-r(1)).^2/(2*r_sig^2));


function replot_quad(quad1,quad2,quad3,quad4,center1,center2,center3,center4,X1,X2,X3,X4,Y1,Y2,Y3,Y4)

c2=center2-center1; c3=center3-center1; c4=center4-center1; 
xmin=min(X1(1,:)); xmax=max(X4(1,:))+max(c3(1),c4(1)); xstep=X1(1,2)-X1(1,1);
ymin=min(Y1(:,1)); ymax=max(Y4(:,1))+max(c3(2),c4(2)); ystep=Y1(2,1)-Y1(1,1);
X=xmin:xstep:xmax; Y=ymin:ystep:ymax;

% prepare zeros
figure(8); hold off;
image(X,Y,zeros(length(X),length(Y))); hold on;


myscale=10;

image(X1(1,:),Y1(:,1),quad1/myscale); hold on;
image(X2(1,:)+c2(1),Y2(:,1)+c2(2),quad2/myscale);
image(X3(1,:)+c3(1),Y3(:,1)+c3(2),quad3/myscale);
image(X4(1,:)+c4(1),Y4(:,1)+c4(2),quad4/myscale);

xlim([0 1800]); ylim([0 1800]);

X;
