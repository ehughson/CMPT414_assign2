I = imread('animals.jpg');
template_img = imread('template_bear.png');
%image_original = img; 

img = rgb2gray(I);
template_img = rgb2gray(template_img);

%First we convert the template image into an edge image using canny edge
%detector

temp2 = edge(template_img,'canny');
figure;
imshow(temp2);
%then we pick a reference point

refx = round(size(temp2,1)/2);
refy = round(size(temp2,2)/2);

%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(temp2), filter1); 
dx = imfilter(double(temp2), filter2); 
[p, q] = size(dx);
angle = size(dx); 


for i = 1:p
    for j = 1:q
        angle(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end


[x, y] = size(temp2);
new_temp = [x, y]; 
tally = 0; 
%disp(temp);
%disp(x);
for m = 1:x
    for n = 1:y
        if(temp2(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end

[x,y] = find(temp2 > 0);
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
    r_table(voting(gradient_index), 1) = reference_x; 
    r_table(voting(gradient_index), 2) = reference_y;
    r_table(voting(gradient_index), 3) = alpha;
end



%create accumulator array --------
img2 = medfilt2(img);
img2 = edge(img2,'sobel', 0.04);

figure;
imshow(img2), title('edge map 1');
%compute the gradients rference
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(img2), filter1); 
dx = imfilter(double(img2), filter2); 
[p, q] = size(dx);
angle2 = size(dx); 


for i = 1:p
    for j = 1:q
        angle2(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end

[x, y] = size(img2);
tally = 0; 
%disp(temp);
%disp(x);
for m = 1:x
    for n = 1:y
        if(img2(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
[a,b] = find(img2 > 0);
%create accumulator array 
accArray= zeros(size(img2)); 
[K, R] = size(img2);

for i = 1:tally
    gradient_index = round(angle2(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        % alpha = r_table(gradient_index, j, 3);
        alpha = -30*pi()/180;
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

index_x = 1; 
index_y = 1; 
prev_index_x = 1;
prev_index_y = 1; 

for i = 1: p
    for j = 1: q
        if MAX <= accArray(i, j)
            MAX = accArray(i, j); 
            prev_index_x = index_x; 
            prev_index_y = index_y;
            index_x = i; 
            index_y = j; 
        end
    end
end

disp(MAX); 
%disp(prev_index_x);
disp(index_x);
disp(index_y); 
        
%xdatatemp = accArray(:, [350:360 318:323]);
%disp(xdatatemp);



figure;
accArray3 = mat2gray(accArray);
imshow(accArray3), title('acc array #1');
plot(index_y,index_x,'g+', 'MarkerSize', 2);


%then we pick a reference point

%temp = medfilt2(temp);
temp3 = edge(template_img,'canny');
refx = round(size(temp3,1)/2);
refy = round(size(temp3,2)/2);

%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(temp3), filter1); 
dx = imfilter(double(temp3), filter2); 
[p, q] = size(dx);
angle = size(dx); 


for i = 1:p
    for j = 1:q
        angle(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end


[x, y] = size(temp3);
tally = 0; 
%disp(x);
for m = 1:x
    for n = 1:y
        if(temp3(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end

[x,y] = find(temp3 > 0);
%create r table
voting = zeros(180, 1); 
r_table2 = zeros(180, 3); 

for i = 1:tally
    reference_x = x(i) - refx;
    reference_y = y(i) - refy; 
    gradient_index = angle(x(i), y(i))+90; 
    %disp(gradient_index); 
    voting(gradient_index) = voting(gradient_index) +1; 
    alpha = atan2d(reference_y, reference_x);
    r_table2(voting(gradient_index), 1) =1.4* reference_x; 
    r_table2(voting(gradient_index), 2) =1.4* reference_y;
    r_table2(voting(gradient_index), 3) = alpha;
end



%create accumulator array --------
img3 = medfilt2(img);
img3 = edge(img3,'sobel', 0.04);


%compute the gradients
filter1 = [1;  -1]; 
filter2 = [1  -1];
dy = imfilter(double(img3), filter1); 
dx = imfilter(double(img3), filter2); 
[p, q] = size(dx);
angle2 = size(dx); 


for i = 1:p
    for j = 1:q
        angle2(i, j) = atan2d(dy(i, j), dx(i, j)); 
    end
end

[x, y] = size(img3);
tally = 0; 
%disp(x);
for m = 1:x
    for n = 1:y
        if(img3(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
[a,b] = find(img3 > 0);
%create accumulator array 
accArray2= zeros(size(img3)); 
[K, R] = size(img3);

for i = 1:tally
    gradient_index = round(angle2(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        % alpha = r_table(gradient_index, j, 3);
        alpha = -105*pi()/180;
        sin_val = sin(alpha); 
        cos_val = cos(alpha);
        nX = round(a(i) - ((r_table2(j, 1)*cos_val) - (r_table2(j, 2)*sin_val)));
        if(nX < 0 || nX > K)
            nX = 0; 
        end
        %disp(nX);
        nY = round(b(i) - ((r_table2(j, 1)*sin_val) + (r_table2(j, 2)*cos_val)));
        if(nY < 0 || nY > R)
            nY = 0; 
        end
        if(nY ~= 0 && nX ~= 0)
            accArray2(nX, nY) = accArray2(nX, nY) +1; 
        end
    end
end



%disp(accArray2(1,1));
[p, q] = size(accArray);



MAX = accArray2(1,1);

index_x2 = 1; 
index_y2 = 1; 
prev_index_x = 1;
prev_index_y = 1; 

for i = 1: p
    for j = 1: q
        if MAX <= accArray2(i, j)
            MAX = accArray2(i, j); 
            prev_index_x = index_x2; 
            prev_index_y = index_y2;
            index_x2 = i; 
            index_y2 = j; 
        end
    end
end

%disp(index_x);
%disp(index_y); 

%disp(MAX);
        
%xdatatemp = accArray2(:, [235:245 255:265]);
%disp(xdatatemp);
figure;
imshow(I);
hold on;
plot(index_y,index_x,'g+', 'MarkerSize', 2);
%hold on;
Circle(index_y, index_x, (1.5 *refy) + 4);
%hold on; 
plot(index_y,index_x,'g+', 'MarkerSize', 2);
%hold on;
Circle(index_y2, index_x2, (1.5 *refy) + 4 );



figure;
accArray4 = mat2gray(accArray2);
imshow(accArray4), title('acc array #2');

%figure;
%imshow(img3);

function Circle(centery, centerx,  r)
angle = 0:pi/50:2*pi; 
d_x = r*cos(angle);
d_y = r*sin(angle);
plot(centery+d_y, centerx+d_x, 'r');
end