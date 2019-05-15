function l = get_laplacian(dim, level)
% return discrete Laplacian operator.
% dim: 1,2,3;
% level:
% - 1d 1(3p)/2(5p)/3(7p)/4(9p);
% - 2d 1(5p)/2(9p);
% - 3d 1(7p)/2(19p)/3(27p)

    switch dim
        case 1
            switch level
                case 1
                    l = [1,-2,1];
                case 2
                    l = [-1/12,4/3,-5/2,4/3,-1/12];
                case 3
                    l = [1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90];
                case 4
                    l = [-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560];
                otherwise
                    error('Level is not defined');
            end
        case 2
            switch level
                case 1
                    l = [0,1,0;1,-4,1;0,1,0];
                case 2
                    l = [1,1,1;1,-8,1;1,1,1]/3;
                otherwise
                    error('Order is not defined');
            end
        case 3
            l = zeros(3,3,3);
            switch level
                case 1
                    l(:,1) = [0,0,0;0,1,0;0,0,0];
                    l(:,2) = [0,1,0;1,-6,1;0,1,0];
                    l(:,3) = [0,0,0;0,1,0;0,0,0];
                case 2
                    l(:,1) = [0,1,0;1,2,1;0,1,0]/6;
                    l(:,2) = [1,2,1;2,-24,2;1,2,1]/6;
                    l(:,3) = [0,1,0;1,2,1;0,1,0]/6;
                case 3
                    l(:,1) = [2,3,2;3,6,3;2,3,2]/26;
                    l(:,2) = [3,6,3;6,-88,6;3,6,3]/26;
                    l(:,3) = [2,3,2;3,6,3;2,3,2]/26;
                otherwise
                    error('Order is not defined');
            end   
        otherwise
            error('Dim is not defined');
    end

end