clear all;
close all;
% mex -I./includes pixhawk_sil_connector.cpp
mex -I./includes pixhawk_sil_connector.cpp -lws2_32