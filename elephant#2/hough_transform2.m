I = imread('animals2.jpg');
template_img = imread('template_elephant.png');
%image_original = img; 

img = rgb2gray(I);
template_img = rgb2gray(template_img);

%First we convert the template image into an edge image using canny edge
%detector


template_img = edge(template_img,'canny');
figure; 
imshow(template_img);
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
for m = 1:x
    for n = 1:y
        if(template_img(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
[x,y] = find(template_img > 0);
%create r table

voting = zeros(180, 1); 

r_table = zeros(180, 3);  

for i = 1:tally
    reference_x = x(i) - refx;
    reference_y = y(i) - refy; 
    gradient_index = angle(x(i), y(i))+90; 
    %disp(gradient_index); 
    voting(gradient_index) = voting(gradient_index) +1; 
    alpha = atan2d(reference_y, reference_x);
    r_table(voting(gradient_index), 1) = 1.4*reference_x; 
    r_table(voting(gradient_index), 2) = 1.4*reference_y;
    r_table(voting(gradient_index), 3) = alpha;
end



%create accumulator array --------
img = medfilt2(img);
img = edge(img,'sobel');
figure; 
imshow(img);
%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(img), filter1); 
dx = imfilter(double(img), filter2); 
[p, q] = size(dx); 
angle = size(dx);

for i = 1:p
    for j = 1:q
        angle(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end

[x, y] = size(img);
 
tally = 0; 
%disp(x);
for m = 1:x
    for n = 1:y
        if(img(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
[a,b] = find(img > 0);
%create r table

accArray= zeros(size(img)); 
[K, R] = size(img);

for i = 1:tally
    gradient_index = round(angle(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        % alpha = r_table(gradient_index, j, 3);
        alpha = 330*pi()/180;
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

        

[p, q] = size(accArray);

MAX = accArray(1,1);

figure; 
imshow(accArray);

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

disp(MAX); 

disp(index_x);
disp(index_y);

xdatatemp = accArray(:, [205:215 75:85]);
%disp(xdatatemp);
figure;
imshow(I);
hold on;
plot(index_y,index_x,'g+', 'MarkerSize', 2);
%hold on;
Circle(index_y, index_x, refy + 4);

function Circle(centery, centerx,  r)
angle = 0:pi/50:2*pi; 
d_x = r*cos(angle);
d_y = r*sin(angle);
plot(centery+d_y, centerx+d_x, 'r');
end
