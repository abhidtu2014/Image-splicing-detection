function [ed] = Edge_detection(I)
ed = edge(I,'sobel',0.027);

