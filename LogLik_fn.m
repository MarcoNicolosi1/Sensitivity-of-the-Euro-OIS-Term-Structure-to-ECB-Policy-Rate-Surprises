function LogLik = LogLik_fn(hyper_par, y, t, T, t_J, N_Estr)


    [Z, d, GG, T_mtx, c, HH, a, P] = SetStateSpaceModel_3F_QF(hyper_par, t, T, t_J,N_Estr);


    LogF = 0; SumSquares = 0;     
    [n_obs,~]= size(y);   

    y = y';
    for i = 1:n_obs
	    v = y(:,i) - Z(:,:,i) * a - d(:,i); 		
        F = Z(:,:,i) * P * Z(:,:,i)' + GG;
        K = T_mtx(:,:,i) * P * Z(:,:,i)' / F;
	    a = T_mtx(:,:,i) * a + c(:,i)+ K * v;		    
        P = T_mtx(:,:,i) * P * T_mtx(:,:,i)' + HH(:,:,i) - K * F * K';
        LogF = LogF + log(det(F));   
        SumSquares = SumSquares + v'*(F\v);
    end
    LogLik = 0.5 * (LogF + SumSquares ); 



end
 