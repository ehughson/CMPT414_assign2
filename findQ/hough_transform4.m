I = imread('letters.png');
template_img = imread('template_Q.png');
%image_original = img; 

img = rgb2gray(I);
template_img = rgb2gray(template_img);
%First we convert the template image into an edge image using canny edge
%detector

template_img = edge(template_img,'canny');
figure; 
imshow(template_img);
%if I used sobel here, it wont detect the third straight oragne Q -- for
%report 

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
new_temp = [x, y]; 
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
%reference https://en.wikipedia.org/wiki/Generalised_Hough_transform

voting = zeros(180, 1); 

r_table = zeros(180, 3); 
r_table2 = zeros(180, 3);


for i = 1:tally
    reference_x = refx - x(i);
    reference_y = refy- y(i); 
    gradient_index = angle(x(i), y(i))+90; 
    %disp(gradient_index); 
    voting(gradient_index) = voting(gradient_index) +1; 
    alpha = atan2d(reference_y, reference_x);
    r_table(voting(gradient_index), 1) = reference_x; 
    r_table(voting(gradient_index), 2) = reference_y;
    r_table(voting(gradient_index), 3) = alpha;
    
    %need to make a second r table becyase the distance is calculated
    %different then when we are dealing with the image straight
    r_table2(voting(gradient_index), 1) = x(i) - refx; 
    r_table2(voting(gradient_index), 2) = y(i) - refy;
    r_table2(voting(gradient_index), 3) = alpha;
end



%create accumulator array --------
img = edge(img,'sobel', 0.024);
figure; 
imshow(img);

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
%disp(x);
for m = 1:x
    for n = 1:y
        if(img(m, n) > 0)
            tally = tally +1;         
            
        end
    end
end
[a,b] = find(img > 0);
%create  accumulator array
accArray= zeros(size(img)); 
[K, R] = size(img);

for i = 1:tally
    gradient_index = round(angle2(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        alpha = r_table( j, 3);
       % alpha = 345*pi()/180;
        sin_val = sin(alpha); 
        cos_val = cos(alpha); 
        %rho = delta_x*cos(alpha) + delta_y*sin(alpha); 
        nX = round(a(i) + r_table( j, 1));
        if(nX < 0 || nX > K)
            nX = 0; 
        end
        %disp(nX);
        nY = round(b(i) + r_table( j, 2));
        if(nY < 0 || nY > R)
            nY = 0; 
        end
        if(nY ~= 0 && nX ~= 0)
            accArray(nX, nY) = accArray(nX, nY) +1; 
        end
    end
end


accArray2= zeros(size(img)); 
for i = 1:tally
    gradient_index = round(angle2(a(i), b(i))+90); 
    %disp(gradient_index); 
    for j = 1:voting(gradient_index)
        alpha = 25*pi()/180;
        sin_val = sin(alpha); 
        cos_val = cos(alpha); 
        nX = round(a(i) - ((r_table2( j, 1)*cos_val) - (r_table2( j, 2)*sin_val)));
        if(nX < 0 || nX > K)
            nX = 0; 
        end
        %disp(nX);
        nY = round(b(i) - ((r_table2( j, 1)*sin_val) + (r_table2(j, 2)*cos_val)));
        if(nY < 0 || nY > R)
            nY = 0; 
        end
         if(nY ~= 0 && nX ~= 0)
            accArray2(nX, nY) = accArray2(nX, nY) +1; 
        end
    end
end



disp(size(accArray));

[p, q] = size(accArray);
[p2, q2] = size(accArray);

MAX = accArray(1,1);

index_x = 1; 
index_y = 1; 

for i = 1: p
    for j = 1: q
        if MAX <= accArray(i, j)
            MAX = accArray(i, j); 
            %disp(accArray(i,j));
            index_x = i; 
            index_y = j; 
        end
    end
end

%disp(MAX);
%disp(index_x); 
%disp(index_y); 
%xdatatemp = accArray(:, [170:180 310:320]);
%disp(xdatatemp);

MAX3 = accArray(1,1);

index_x3 = 1; 
index_y3 = 1; 
for i = 1: p
    for j = 1: q
        if (MAX3 <= accArray(i, j) && i ~= index_x && j ~= index_y)
            MAX3 = accArray(i, j); 
            %disp(accArray(i,j));
            index_x3 = i; 
            index_y3 = j; 
        end
    end
end

%disp(MAX3);
%disp(index_x3); 
%disp(index_y3); 
%xdatatemp = accArray(:, [193:196 570:575]);
%disp(xdatatemp);

MAX2 = accArray2(1,1);

index_x2 = 1; 
index_y2 = 1;
for i = 1: p2
    for j = 1: q2
        if MAX2 <= accArray2(i, j)
            MAX2 = accArray2(i, j); 
            %disp(accArray(i,j));
            index_x2 = i; 
            index_y2 = j; 
        end
    end
end

%disp(MAX2);
%disp(index_x2); 
%disp(index_y2); 
%xdatatemp = accArray2(:, [268:272 185:190]);
%disp(xdatatemp);




MAX4 = accArray(1,1);

index_x4 = 1; 
index_y4 = 1; 
for i = 1: p
    for j = 1: q
        if (MAX4 <= accArray(i, j)&& i ~= index_x && j ~= index_y && i ~= index_x3 && j ~= index_y3 )
            MAX4 = accArray(i, j); 
            %disp(accArray(i,j));
            index_x4 = i; 
            index_y4 = j; 
        end
    end
end

%disp(MAX4);
%disp(index_x4); 
%disp(index_y4); 
%xdatatemp = accArray(:, [302:307 242:247]);
%disp(xdatatemp);






figure;
imshow(I);
hold on;
plot(index_y,index_x,'g+', 'MarkerSize', 2);
hold on;
Circle(index_y3, index_x3, refy + 4);
hold on;
Circle(index_y, index_x, refy + 4);
hold on;
%Circle(366, 260, refPointy, 3);
hold on;
Circle(index_y2, index_x2, refy + 4);
hold on;
Circle(index_y4, index_x4, refy + 4);


figure; 
accArrayresult = mat2gray(accArray);
imshow(accArrayresult) ,title('accArray1');

figure; 
accArrayresult2 = mat2gray(accArray2);
imshow(accArrayresult2) ,title('accArray1');

function Circle(centery, centerx,  r)
angle = 0:pi/50:2*pi; 
d_x = r*cos(angle);
d_y = r*sin(angle);
plot(centery+d_y, centerx+d_x, 'r');
end