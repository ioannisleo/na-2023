function xr = newtonraphson(f,x0,es,fs,s)
% Υπολογισμός ριζών της συνάρτησης f με χρήση της μεθόδου Newton-Raphson.
%   f  : συμβολική συνάρτηση
%   x0 : αρχικά σημεία
%   es : ακρίβεια σύγκλισης 
%   fs : ακρίβεια τιμής
%   s  : αποθήκευση

    MAXITER = 100;
    dim = length(argnames(f));
    J = matlabFunction(inv(jacobian(f)),Vars={argnames(f).'});
    f = matlabFunction(f,Vars={argnames(f).'});
    
    xr = NaN(dim,length(x0));
    for i = 1:length(x0)
        % υπολογισμός για κάθε αρχικό σημείο
        x = NaN(dim,MAXITER);
        fx = NaN(dim,MAXITER);
        e = NaN(dim,MAXITER);      
    
        % αρχικοποίηση
        x(:,1) = x0(:,i);
        e(:,1) = Inf(dim,1);
        fx(:,1) = f(x0(:,i));
        
        n = 1;
        while n < MAXITER && ((any(es < e(:,n))) || (fs < norm(fx(:,n))))
            % επαναληπτική προσέγγιση
            x(:,n+1) = x(:,n) - J(x(:,n))*fx(:,n);
            e(:,n+1) = abs(x(:,n+1) - x(:,n));
            fx(:,n+1) = f(x(:,n+1));
            n = n + 1;
        end
        
        if n < MAXITER
            % επιτυχής εύρεση
            xr(:,i) = x(:,n);

            if s == 1
                % αποθήκευση
                results = table(x(:,1:n)',fx(:,1:n)', ...
                    e(:,1:n)',VariableNames={'x','f','e'});
                save(sprintf('results%d.mat',i), "results");
            end
        end

    end
end
