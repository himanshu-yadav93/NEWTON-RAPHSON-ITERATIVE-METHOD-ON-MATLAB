% Program for Newton-Raphson Load Flow Analysis..
% Ref Book: Power System Analysis by Hadi Saadat..

busdata = busdata6();       % Calling "busdata6.m" for bus data.
Y = ybusppg();              % Calling ybusppg.m to get Bus Admittance Matrix..
BMva = 100;                 % Base MVA..
bus = busdata(:,1);            % Bus Number..
type = busdata(:,2);           % Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busdata(:,3);              % Specified Voltage..
del = zeros(length(V),1);   % Voltage Angle..
Pg = busdata(:,5)/BMva;        % PGi..
Qg = busdata(:,6)/BMva;        % QGi..
Pl = busdata(:,7)/BMva;        % PLi..
Ql = busdata(:,8)/BMva;        % QLi..
Qmin = busdata(:,9)/BMva;      % Minimum Reactive Power Limit..
Qmax = busdata(:,10)/BMva;     % Maximum Reactive Power Limit..
nbus = max(bus);            % To get no. of buses..
P = Pg - Pl;                % Pi = PGi - PLi..
Q = Qg - Ql;                % Qi = QGi - QLi..
Psp = P;                    % P Specified       
Qsp = Q;                    % Q Specified
G = real(Y);                % Conductance..
B = imag(Y);                % Susceptance..

pv = find(type == 2 | type == 1);       % Index of PV Buses..
pq = find(type == 3);                   % Index of PQ Buses..
npv = length(pv);                       % Number of PV buses..
npq = length(pq);                       % Number of PQ buses..

Tol = 1;  
Iter = 1;           % iteration starting
while (Tol > 1e-5)   % Iteration starting..
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculate P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    % Checking Q-limit violations..
    if Iter <= 7 && Iter > 2    % Only checked up to 7th iterations..
        for n = 2:nbus
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];         % Jacobian Matrix
    
    X=inv(J)*M;                 % INV(J) x M, Correction Vector..
    dTh = X(1:nbus-1);          % Change in Voltage Angle..
    dV = X(nbus:end);           % Change in Voltage Magnitude..
    
    % Update State Vectors (Voltage Angle & Magnitude)
    del(2:nbus) = dTh + del(2:nbus);
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);
            k = k+1;
        end
    end
    Iter = Iter + 1;
    Tol = max(abs(M));
end
Del = 180/pi*del;   % Convert radians to degrees

 %% Load Flow Solution 
disp('----------------------------------------');
fprintf( 'Number Of Ieration %3g \n',Iter)

disp('----------------------------------------');
disp('  Newton Raphson Loadflow Solution    ');
disp('----------------------------------------');
disp(' |Bus |   |Voltage|    |Angle |');
disp(' | No.|   |pu     |    |Degree|');
disp('----------------------------------------');
for m=1:nbus
    fprintf(' %3g   ' ,m);
    fprintf(' %8.4f    ' ,V(m));
    fprintf(' %8.4f  ' ,Del(m));
    fprintf('\n');
end

 