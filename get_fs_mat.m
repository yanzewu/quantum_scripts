
% Return a derivation matrix using finite stencil
function ret = get_fs_mat(n, order, acc, ret_mat)

    if ~exist('ret_mat', 'var')
        ret_mat = true;
    end
    
    if n == 0
        ret_mat = false;
    end

    switch order
        case 1
            filler_mat = [-1/2,0,1/2,0,0,0,0;
                            1/12,-2/3,0,2/3,-1/12,0,0;
                            -1/60,3/20,-3/4,0,3/4,-3/20,1/60];
        case 2
            filler_mat = [1,-2,1,0,0,0,0,0,0;
                            -1/12,4/3,-5/2,4/3,-1/12,0,0,0,0;
                            1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90,0,0;
                            -1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560];
        case 4
            filler_mat = [1,-4,6,-4,1,0,0,0,0,0,0;
                        -1/6,2,-13/2,28/3,-13/2,2,-1/6,0,0,0,0;
                        7/240,-2/5,169/60,-122/15,91/8,-122/15,169/60,-2/5,7/240,0,0;
                        -82/15120,1261/15120,-9738/15120,52428/15120,-140196/15120,192654/15120,-140196/15120,52428/15120,-9738/15120,1261/15120,-82/15120];

        otherwise
            error('Order is not defined');
    end

    % Load stencil vector from storage;
    filler = filler_mat(acc/2,1:order+acc-1+mod(order,2));
    filler = filler(1:order+acc-1+mod(order,2));

    if ret_mat % return a matrix

        half_sz_filler = (length(filler)-1)/2;

        fs_mat = zeros(n,n);
        for j = half_sz_filler+1:n-half_sz_filler
            fs_mat(j,j-half_sz_filler:j+half_sz_filler) = filler;
        end
        
        for j = 1:half_sz_filler
            fs_mat(j,1:j+half_sz_filler) = filler(half_sz_filler-j+2:length(filler));
            fs_mat(n-j+1,n-j+1-half_sz_filler:n) = filler(1:j+half_sz_filler);
        end

        ret = fs_mat;
    else
        ret = filler;
    end
end