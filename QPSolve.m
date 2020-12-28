classdef QPSolve
	properties
		xdim
        udim
        dyn
        clf
        cbf_list
        u_lim
        p1_cost
        p2_cost
        verbose
        u_cost
        u_prev_cost
        K
        ksig 
        max_var
        mu_qp_prev
        P
        A
        A0
        G
        res
        max_error
	end
	methods
		function obj = QPSolve(dyn,cbf_list,clf,u_lim,u_cost,u_prev_cost,p1_cost,p2_cost,verbose)
			obj.xdim = dyn.xdim;
			obj.udim = dyn.udim;
			obj.dyn = dyn;
			obj.cbf_list = cbf_list;
			obj.clf = clf;
			obj.u_lim = u_lim;
			obj.p1_cost = p1_cost;
			obj.p2_cost = p2_cost;
			obj.verbose = verbose;
			obj.u_cost = u_cost;
			obj.u_prev_cost = u_prev_cost;
			obj.K = 0.0;
			obj.ksig = 1.0;
			obj.max_var = 1.0;
			obj.mu_qp_prev = zeros(floor(obj.xdim/2),1,'double');
			obj.P = eye(obj.xdim,'double');
			obj.A = zeros(obj.xdim,obj.xdim,'double');
			obj.A0 = double([zeros(floor(obj.xdim/2)),eye(floor(obj.xdim/2)); ...
							   zeros(floor(obj.xdim/2)),zeros(floor(obj.xdim/2))]);
			obj.G = double([zeros(floor(obj.xdim/2));eye(floor(obj.xdim/2))]);
			obj.res = NaN;
			obj.max_error = 1.0;
		end
	
		function update_ricatti(obj,A)
			obj.A = A;
			Q = eye(obj.xdim,'double');
			[X,K,L] = icare(A,B,Q,R,S,E,G)
			[obj.P,K,L] = icare(obj.A,)
			obj.P = sl.solve_continuous_are(obj.A,np.zeros((obj.xdim,obj.xdim), dtype=np.float32),Q,np.eye(obj.xdim, dtype=np.float32))
		end
		% untest
		function result = solve(obj,x,x_d,mu_d,sigDelta)
			sigDelta = sigDelta * obj.ksig;
			sigDelta = min(max(sigDelta,0.0),obj.max_var); % np.clip(sigDelta,0.0,obj.max_var)

			% build Q and p matrices to specify minimization expression
			Q = diag(np.append(np.append(double(ones(floor(obj.xdim/2)))*(obj.u_cost + obj.u_prev_cost),obj.p1_cost),obj.p2_cost))
			obj.Q = sparse(Q);
			obj.p = 2*[-obj.mu_qp_prev*obj.u_prev_cost,0,0];

			% #error dynamics for clf
			e = x(1:length(x)-1,:)-x_d(1:length(x)-1,:); %x[:-1,:]-x_d[:-1,:]
			e = e(e>obj.max_error)=obj.max_error;% np.clip(e,-obj.max_error,obj.max_error)
			e = e(e<-obj.max_error)=-obj.max_error;
			eTPG = e.T*(obj.P*obj.G);
			
			G_dyn = [[eTPG,1,0]]; %#append 1 for clf < d
			Gsig = obj.G*sigDelta;
			GssG = Gsig*Gsig.T;
			obj.trGssGP = trace(GssG*obj.P);
			h_dyn = -1 * ( -0.5*np.matmul(e.T*(Q*e)) ...
						+ 0.5*(e.T*(obj.P*e)) / obj.clf.epsilon ...
						+ 0.5*obj.trGssGP);

			% # build constraints for barriers
			N_cbf = length(obj.cbf_list);
			G_cbf = double(zeros(N_cbf,floor(obj.xdim/2)+2));
			h_cbf = double(zeros(N_cbf,1);
			A0x_Gmud = obj.A0*x(1:length(x)-1,:)) + obj.G*mu_d;
			GssG_22 = GssG(3:length(x),3:length(x));
			for i = 1:length(obj.cbf_list)
				cbf = obj.cbf_list(i);
				[h_x, dB, d2B] = cbf.get_B_derivatives(x);
				G_cbf(i,:) = [einsum('ik,kj-> ij',dB,obj.G),0,1];
				trGssGd2B = trace(einsum('ik,kj-> ij',GssG_22,d2B));
				h_cbf(i,:) = -1 * (einsum('ik,kj-> ij',dB,A0x_Gmud) ...
									- cbf.gamma * h_x ...
									+ 0.5*trGssGd2B);
			end
			% # build constraints for control limits
			ginv = inv(obj.dyn.g(x));
			l_ctrl = ginv*( mu_d - obj.dyn.f(x));
			A_ctrl = ginv;

			G_ctrl = double(zeros(obj.udim*2,floor(obj.xdim/2)+2));
			h_ctrl = double(zeros(obj.udim*2,1));
			for i = 1:length(obj.udim)
				G_ctrl(i*2,1:floor(obj.xdim/2)+1] = - A_ctrl(i,:);
				h_ctrl(i*2) = - obj.u_lim(i,1) + l_ctrl(i);
				G_ctrl(i*2+1,1:floor(obj.xdim/2)+1] = A_ctrl(i,:);
				h_ctrl(i*2+1) = obj.u_lim(i,2) - l_ctrl(i);
			end
			% # stack into one matrix and vector
			G = cat(1,G_dyn,G_cbf,G_ctrl);
			h = cat(1,h_dyn,h_cbf,h_ctrl);

			obj.G_csc = sparse(G);
			obj.h = h;

			% # dummy lower bound
			l = double(ones(h.shape))*Inf * -1;

			% #Solve QP
			obj.prob = osqp;
			exception_called = false;
			mu_bar = double(zeros(obj.xdim+1));
			% obj.prob.setup(P=obj.Q, q=obj.p, A=obj.G_csc, l=l, u=obj.h, verbose=obj.verbose)
			obj.prob.setup(obj.Q, obj.p, obj.G_csc, l, obj.h, obj.verbose);
			obj.res = obj.prob.solve();
			mu_bar = obj.res.x;
			if any(isnan(mu_bar))
				mu_bar = zeros(obj.xdim+1);
				obj.res = nan;
				disp("QP failed!");
			end

			obj.mu_qp_prev = [mu_bar[1:floor(obj.xdim/2)+1)]';

			obj.V =obj.clf.V(x,x_d);

			% if obj.verbose:
				% print('z_ref: ', x.T)
				% print('z_des: ', x_d.T)
				% print('u_lim', obj.u_lim)
				% print('V: ', obj.V)
				% print('Q:', Q)
				% print('p:', np.array(obj.p))
				% print('G_dyn:', G_dyn)
				% print('h_dyn:', h_dyn)
				% print('trGssGP',obj.trGssGP)
				% if h_cbf.shape[0] < 10:
					% print('G_cbf:', G_cbf)
					% print('h_cbf:', h_cbf)
				% print('G_ctrl:', G_ctrl)
				% print('h_ctrl:', h_ctrl)
				% print('result:', mu_bar)
				% print('G_all:',obj.G_csc)
				% print('h_all:',obj.h)

			return obj.mu_qp_prev
		end
	end
end	