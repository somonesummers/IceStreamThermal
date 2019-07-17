clear all
clc

for i = 1:4
    cases = ['dnB';'upB';'ISE';'profileI'];
    cases2 = ['dnWhillans';'upWhillans';'macAyeal';'bindschadler'];
    offset = [0;0;0;0];
    load digitizedVelocities.mat
    hold on
    evaluation = sprintf('plot((-%s(:,1)+%d)/1e3,%s(:,2),''.'')',cases(i,:),offset(i,:),cases(i,:),i);
    eval(evaluation)
end