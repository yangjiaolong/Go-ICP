% Simple demos for the GoICP code
% Demo1 uses 100 points randomly drawn in [-1,1]^3
% Demo2 uses the Stanford bunny data normalized in [-1,1]^3
% Gaussian noise is added to the all points
%
% Jiaolong Yang <yangjiaolong@gmail.com>
% Feb 13, 2014


% data = 'random points';
data = 'bunny';

if ispc()
    cmd = 'GoICP_vc2012.exe';
else
    cmd = './GoICP';
end

if strcmp(data, 'random points')
    cmd = [cmd ' model_rand.txt data_rand.txt 100 config.txt output.txt'];
    model = readpoints('model_rand.txt');
    data = readpoints('data_rand.txt');
else
    cmd = [cmd ' model_bunny.txt data_bunny.txt 1000 config.txt output.txt'];
    model = readpoints('model_bunny.txt');
    data = readpoints('data_bunny.txt');
end

system(cmd);

file = fopen('output.txt', 'r');
t = fscanf(file, '%f', 1);
R = fscanf(file, '%f', [3,3])';
T = fscanf(file, '%f', [3,1]);
fclose(file);

figure;
subplot(1,2,1);
plot3(model(1,:), model(2,:), model(3,:), '.r');
hold on;
plot3(data(1,:),  data(2,:),  data(3,:), '.b');
hold off; axis equal; title('Initial Pose');
subplot(1,2,2);
data_ = bsxfun(@plus, R*data, T);
plot3(model(1,:), model(2,:), model(3,:), '.r');
hold on;
plot3(data_(1,:), data_(2,:), data_(3,:), '.b');
hold off; axis equal;  title('Result');
