fid = fopen('latticeMatch_gr.csv', 'wt');
a = [4.34052 0];
b = [4.34052*cos(deg2rad(60)) 4.34052*sin(deg2rad(60))];

matchSize = 14.75707;

syms b_scaleX b_scaleY;

fprintf(fid,"Lattice const,a1,b1,a2,b2,theta,match\n");

for i = -10:10
    for j = -10:10
        for k = -10:10
            for l = -10:10
                if ~(i == 0 && j == 0 && k == 0 && l == 0)
                    a_prime = i.*a + j.*b;
                    b_prime = k.*a + l.*b;
                    % if(norm(a_prime) <= 16)
                    % disp(strcat("norm_a ",string(norm(a_prime)),", norm_b ",string(norm(b_prime))));
                    % end
                    if(norm(a_prime) <= matchSize+3 && round(norm(a_prime),3) == round(norm(b_prime),3))
                        theta = acosd(dot(a_prime,b_prime)/norm(a_prime)^2);
                        % disp(strcat("n_m, theta: ", string(theta)));
                        if(round(theta) == 60)
                            disp(strcat(string(i), "*a +", string(j), "*b, ", string(k), "*a +", string(l), "*b, mismatch: ", string(norm(a_prime)/matchSize)));
                            fprintf(fid,strcat(string(norm(a_prime)),",",string(i),",", string(j),",", string(k),",", string(l),",",string(theta),",", string(norm(a_prime)/matchSize),"\n"));
                        end
                    end
                end
            end
        end
    end
end

fclose('all');