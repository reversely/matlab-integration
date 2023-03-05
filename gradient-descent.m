f = @(x) sin(x(1))*cos(1.5*x(2)) - cos(x(1))*sin(x(2)) - sin(1.1*x(1))*sin(x(2));
x1 = [0,-1]; h = 0.1; alpha = 0.2; tol = 1e-5; max_iter = 100;
X = gradient_descent(f,x1,alpha,h,tol,max_iter);
F = @(x,y) sin(x).*cos(1.5*y) - cos(x).*sin(y) - sin(1.1*x).*sin(y);
points = linspace(-pi,pi,50);
[x,y] = meshgrid(points,points);
z = F(x,y);
contour(x,y,z), hold on
for k=1:length(X)
    plot(X(k,1),X(k,2),'r.')
end


%% Function for approximation of gradient vector

function df = gradient(f,a,h)
n = length(a);
df = zeros(1,n);
for i=1:n
    hi = zeros(1,n);
    hi(i) = h;
    df(i) = (f(a + hi) - f(a - hi))/(2*h);
end
end

function X = gradient_descent(f,x1,alpha,h,tol,max_iter)
tstart=tic;
n = length(x1);
X = zeros(max_iter+1,n);
X(1,:) = x1;

for k=1:max_iter
    grad = gradient(f,X(k,:),h);
    X(k+1,:) = X(k,:) - alpha*grad;
    grad = gradient(f,X(k+1,:),h);
    a = X(k+1,:);
    ngrad = norm(grad);

if (ngrad < tol) || (k==max_iter)

        % calculating Hessian
        H = zeros(n,n);
        for i=1:n
            for j=(i+1):n
                hi = zeros(1,n); hi(i) = h;
                hj = zeros(1,n); hj(j) = h;
                H(i,j) = (f(a + hi + hj) - f(a + hi - hj) - f(a - hi + hj) + f(a - hi - hj))/(4*h^2);
                H(j,i) = H(i,j);
            end
        end
        for i=1:n
            hi = zeros(1,n); hi(i) = h;
            H(i,i) = (f(a + hi) - 2*f(a) + f(a - hi))/(h^2);
        end

        % compute second derivative test
        evals = eig(H);
        telapsed = toc(tstart);

        if all(evals>0)
            disp(['Found local minimum in ', num2str(telapsed), ' seconds and ', num2str(k), ' iterations.'])
            break

        elseif all(evals<0)
            disp(['Found local maximum in ', num2str(telapsed), ' seconds and ', num2str(k), ' iterations.'])
            break

        elseif any(evals==0)
            disp(['Found local saddlepoint (or degenerate point) in ', num2str(telapsed), ' seconds and ', num2str(k), ' iterations.'])
            X = "NaN"
            break

        elseif k==max_iter
            X = "NaN"
            disp("Exceeded maximum iterations. No solution found.")

        else
            disp(["Found local saddlepoint (or degenerate point) in ", num2str(telapsed), " seconds and ", num2str(k), " iterations."])
            X = "NaN"
            break


    end

end
end

end
