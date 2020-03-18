I = imread('animals.jpg');
template_img = imread('template_elephant.png');
%image_original = img; 

img = rgb2gray(I);
template_img = rgb2gray(template_img);

%First we convert the template image into an edge image using canny edge
%detector

template_img = edge(template_img,'sobel');

figure; 
imshow(template_img), title('edge map for sobel'); 

%then we pick a reference point

refx = round(size(template_img,1)/2);
refy = round(size(template_img,2)/2);

%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(template_img), filter1); 
dx = imfilter(double(template_img), filter2); 
[p, q] = size(dx);
angle = size(dx); 


for i = 1:p
    for j = 1:q
        angle(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end

[x, y] = size(template_img);
 
tally = 0; 
%disp(x);

%take tally of all edge points
for m = 1:x
    for n = 1:y
        if(template_img(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
%disp(size(tally));
%disp(tally);
[x,y] = find(template_img > 0);
%disp(size(x));


%create r table

voting = zeros(180, 1);
r_table = zeros(180, 3); 


for i = 1:tally
    reference_x = x(i) - refx;
    reference_y = y(i) - refy; 
    gradient_index = round(angle(x(i), y(i)))+90; 
    %disp(gradient_index); 
    voting(gradient_index) = voting(gradient_index)+1; 
    %disp(voting(gradient_index)); 
    %compute alpha which is the angle between the horizontal direction and
    %the reference point
    alpha = atan2d(reference_y, reference_x);
    r_table(voting(gradient_index), 1) = reference_x; 
    r_table(voting(gradient_index), 2) = reference_y;
    r_table(voting(gradient_index), 3) = alpha;
end



%create accumulator array --------
%img = medfilt2(img);
img = edge(img,'canny', 0.45);

figure; 
imshow(img), title('edge map'); 


%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(img), filter1); 
dx = imfilter(double(img), filter2); 
[p, q] = size(dx);
angle2 = size(dx); 


for i = 1:p
    for j = 1:q
        angle2(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end

[x, y] = size(img);
new_img = [x, y]; 
tally = 0; 
%disp(temp);
%disp(x);
for m = 1:x
    for n = 1:y
        if(img(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end

[a,b] = find(img > 0);
%create accumulator array
disp(size(img)); 
accArray= zeros(size(img)); 
[K, R] = size(img);

for i = 1:tally
    gradient_index = round(angle2(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        % alpha = r_table(gradient_index, j, 3);
        alpha = pi()/2;
        sin_val = sin(alpha); 
        cos_val = cos(alpha);
        nX = round(a(i) - ((r_table( j, 1)*cos_val) - (r_table( j, 2)*sin_val)));
        if(nX < 0 || nX > K)
            nX = 0; 
        end
        %disp(nX);
        nY = round(b(i) - ((r_table( j, 1)*sin_val) + (r_table( j, 2)*cos_val)));
        if(nY < 0 || nY > R)
            nY = 0; 
        end
        if(nY ~= 0 && nX ~= 0)
            accArray(nX, nY) = accArray(nX, nY) +1; 
        end
    end
end

        
figure; 
accArray2 = mat2gray(accArray);
imshow(accArray2); 

%disp(accArray);
xdatatemp = accArray(:, [220:250 400:420]);
%disp(xdatatemp);
[p, q] = size(accArray);
%disp(size(accArray)); 

%find peak
MAX = accArray(1,1);

index_x = 1; 
index_y = 1; 

for i = 1: p
    for j = 1: q
        if MAX <= accArray(i, j)
            MAX = accArray(i, j); 
            index_x = i; 
            index_y = j; 
        end
    end
end

disp(index_x); 
disp(index_y); 

figure;
imshow(I);
hold on;
plot(index_y,index_x,'g+', 'MarkerSize', 2);
hold on;
Circle(index_y, index_x, refy + 4);

function Circle(centery, centerx,  r)
angle = 0:pi/50:2*pi; 
d_x = r*cos(angle);
d_y = r*sin(angle);
plot(centery+d_y, centerx+d_x, 'r');
end








