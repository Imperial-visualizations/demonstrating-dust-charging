/*jshint esversion: 6 */
$(window).on('load', function() {//main 
    function find_W_H(a_0,z,root_prec){

        function f_x(x_n,z){
            return x_n*Math.exp(x_n) - z;
        }
        function f_x_first_derv(x_n){
            return Math.exp(x_n)*(1 + x_n);
        }
        function f_x_second_derv(x_n){
            return Math.exp(x_n)*(2 + x_n);
        }
        function calc_x_plus(x_n,z){
            let f_x_0 = f_x(x_n,z);
            let f_x_1 = f_x_first_derv(x_n);
            let f_x_2 = f_x_second_derv(x_n);
            let x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(f_x_1**2)-f_x_0*f_x_2));
            return x_plus;
        }

        let a_n = a_0;
        let a_plus = calc_x_plus(a_0,z);

        while(Math.abs(a_plus - a_n) > root_prec){
            //console.log(a_n);
            a_n = a_plus;
            a_plus = calc_x_plus(a_n,z);
            //console.log(a_plus); 
        }
        W_0 = (a_n + a_plus)/2;
        return W_0; 
    }
    function erf(x){

        let a_1 = 0.254829592;
        let a_2 = -0.284496736;
        let a_3 = 1.421413741;
        let a_4 = -1.453152027;
        let a_5 = 1.061405429;
        let p = 0.3275911;
        let t = 1/(1 + p*x);
        let val = 1 - (a_1*t + a_2*t**2 + a_3*t**3 + a_4*t**4 + a_5*t**5) * Math.exp(-1*(x**2));

        return val;
    }

    function ABR_calc(){
        function find_phi(a_0,J,gamma,root_prec){
            function f_x(x_n,J,gamma){
                return (4*(x_n**1.5)*(2*x_n-3)*(2*x_n+1))/((2*x_n-1)**3) - J/gamma;
            }
            function f_x_first_derv(x_n){
                return (x_n**3.5-2.5*x_n**2.5 + 4.75*x_n**1.5+1.125*x_n**0.5)/((0.5-x_n)**4);
            }
            function f_x_second_derv(x_n){
                return (-0.5*x_n**4+2*x_n**3-8.75*x_n**2-7.5*x_n-0.28125)/(((x_n-0.5)**5)*(x_n**0.5));
            }
            function calc_x_plus(x_n,J,gamma){
                let f_x_0 = f_x(x_n,J,gamma);
                let f_x_1 = f_x_first_derv(x_n);
                let f_x_2 = f_x_second_derv(x_n);
                let x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(f_x_1**2)-f_x_0*f_x_2));

                return x_plus;
            }
        
            let a_n = a_0;
            let a_plus = calc_x_plus(a_0,J,gamma);
        
            while(Math.abs(a_plus - a_n) > root_prec){
                a_n = a_plus;
                a_plus = calc_x_plus(a_n,J,gamma);
            }

            phi_b = (a_n + a_plus)/2;
            return phi_b; 
        }

        function boundry_u(a_0,J,g,root_prec){
            let phi_calc = find_phi(a_0,J,g,root_prec);
            return phi_calc;
        }

        function boundry_rho(phi_b,J){
            let rho_b = ((J**(1/2)*math.exp(phi_b/2))/(phi_b**(1/4)));
            return rho_b;
        }

        function boundry_v(init_rho,init_phi,J){
            let v = (2*init_rho*init_phi**(3/2)*math.exp(-init_phi))/(J*(init_phi-1/2));
            return v;
        }
        function inital_conditions(a_0,j,g,root_prec){
            let inital_u = boundry_u(a_0,j,g,root_prec);
            let inital_rho = boundry_rho(inital_u,j);
            let inital_v = boundry_v(inital_rho,inital_u,j);
            return [inital_rho, inital_u, inital_v];
        } 
        function RK45(rho,u,v,h,J,eps){

            let a_2 = 1/5
            let a_3 = 3/10
            let a_4 = 3/5
            let a_5 = 1
            let a_6 = 7/8
    
            let c_1 = 37/378
            let c_2 = 0
            let c_3 = 250/621
            let c_4 = 125/594
            let c_5 = 0
            let c_6 = 512/1771
            
            let c_1_s = 2825/27648
            let c_2_s = 0
            let c_3_s = 18575/48384
            let c_4_s = 13525/55296
            let c_5_s = 277/14336
            let c_6_s = 1/4
            
            let b_21 = 1/5
            let b_31 = 3/40
            let b_32 = 9/40
            let b_41 = 3/10
            let b_42 = -9/10
            let b_43 = 6/5
            let b_51 = -11/54
            let b_52 = 5/2
            let b_53 = -70/27
            let b_54 = 35/27
            let b_61 = 1631/55296
            let b_62 = 175/512
            let b_63 = 575/13824
            let b_64 = 44275/110592
            let b_65 = 253/4096
    
            function f_der_45(rho,beta,J){
                let f = [,];
                let u = beta[0];
                let v = beta[1];
                f[0] = v;
                f[1] = rho**(-2) * J* u**(-1/2)   -   math.exp(-u)  -   2*v*rho**(-1);
                return f;
            }
            function k_1_45(rho,beta,h,J){
                let k_1 = math.multiply(h,f_der_45(rho,beta,J));
                return k_1;
            }
            function k_2_45(rho,beta,h,J,k_1){
                beta = math.add(beta,math.multiply(b_21,k_1));
                rho = rho+a_2*h; 
                let k_2 = math.multiply(h,f_der_45(rho,beta,J));
                return k_2;
            }
            function k_3_45(rho,beta,h,J,k_1,k_2){
                beta = math.add(beta,math.multiply(b_31,k_1),math.multiply(b_32,k_2));
                rho = rho+a_3*h; 
                let k_3 = math.multiply(h,f_der_45(rho,beta,J));
                return k_3;
            }
            function k_4_45(rho,beta,h,J,k_1,k_2,k_3){
                beta = math.add(beta,math.multiply(b_41,k_1),math.multiply(b_42,k_2),math.multiply(b_43,k_3));
                rho = rho+a_4*h; 
                let k_4 = math.multiply(h,f_der_45(rho,beta,J));
                return k_4;
            }
            function k_5_45(rho,beta,h,J,k_1,k_2,k_3,k_4){
                beta = math.add(beta,math.multiply(b_51,k_1),math.multiply(b_52,k_2),math.multiply(b_53,k_3),math.multiply(b_54,k_4));
                rho = rho+a_5*h; 
                let k_5 = math.multiply(h,f_der_45(rho,beta,J));
                return k_5;
            }
            function k_6_45(rho,beta,h,J,k_1,k_2,k_3,k_4,k_5){
                beta = math.add(beta,math.multiply(b_61,k_1),math.multiply(b_62,k_2),math.multiply(b_63,k_3),math.multiply(b_64,k_4),math.multiply(b_65,k_5));
                rho = rho+a_6*h; 
                let k_6 = math.multiply(h,f_der_45(rho,beta,J));
                return k_6;
            }
    
            function RK_fomulae(rho,u,v,h,J){
                let beta = [u,v];
                let k1 = k_1_45(rho,beta,h,J);
                let k2 = k_2_45(rho,beta,h,J,k1);
                let k3 = k_3_45(rho,beta,h,J,k1,k2);
                let k4 = k_4_45(rho,beta,h,J,k1,k2,k3);
                let k5 = k_5_45(rho,beta,h,J,k1,k2,k3,k4);
                let k6 = k_6_45(rho,beta,h,J,k1,k2,k3,k4,k5);
                let beta_5 = math.add(beta,math.multiply(c_1,k1),math.multiply(c_2,k2),math.multiply(c_3,k3),math.multiply(c_4,k4),math.multiply(c_5,k5),math.multiply(c_6,k6));
                let beta_4 = math.add(beta,math.multiply(c_1_s,k1),math.multiply(c_2_s,k2),math.multiply(c_3_s,k3),math.multiply(c_4_s,k4),math.multiply(c_5_s,k5),math.multiply(c_6_s,k6));
                return [beta_5,beta_4]
            }
    
            function step(rho,u,v,h,J,eps){
                let delta_0 = math.abs(math.multiply(eps,[u,v]));
                let S = 0.95;
                let h_new;
                let rk = RK_fomulae(rho,u,v,h,J)
                let beta_5 = rk[0];
                let beta_4 = rk[1];
                let delta_1 = math.abs(math.subtract(beta_5,beta_4));
                let R = Math.max(...math.abs(math.dotDivide(delta_1,delta_0)));
    
                if(R > 1){
                    let h_loop = math.multiply(S*h,Math.pow(R,-0.25));
                    while(R > 1){
                        rk = RK_fomulae(rho,u,v,h_loop,J);
                        beta_5 = rk[0];
                        beta_4 = rk[1];
                        delta_1 = math.subtract(beta_5,beta_4);
                        R = Math.max(...math.abs(math.dotDivide(delta_1,delta_0)));
                    }
                    h_new = h_loop;
                }else{
                    h_new = math.multiply(S*h,Math.pow(R,-0.2));
                }
                let new_rho = rho+h_new;
                let new_cond = [new_rho,beta_5[0],beta_5[1],h_new];
                return new_cond;
            }
            let new_cond = step(rho,u,v,h,J,eps);
            return new_cond;
        }

        function results_eta(a_0,n,j,g,root_prec,eps){

            let current_step = inital_conditions(a_0,j,g,root_prec);
            let h = -1*current_step[0]/n;
            current_step.push(h);
            let phi_data = [];
            let rho_data = [];
            let h_data = [];
            
            rho_data.push(current_step[0]);
            phi_data.push(current_step[1]);
            h_data.push(current_step[3]);

            while(current_step[0] > Math.abs(h)){
                current_step = RK45(current_step[0],current_step[1],current_step[2],current_step[3],j,eps);//rho,u,v,h,J,eps
                rho_data.push(current_step[0]);
                phi_data.push(current_step[1]);
                h_data.push(current_step[3]);
            }

            return [rho_data,phi_data,h_data];
        }

        function produce_plot_floating(){
            let g = 1e3;
            let rootprec = 10**(-12);
            let a_0 = 0.1;
            let n = 5000

            let eps = 10**parseFloat(document.getElementById('EpsController').value); 
            let J = 10**parseFloat(document.getElementById('JController').value); 
            let Z = parseFloat(document.getElementById('ZController').value); 
            let alpha = (Z**0.5)*(1836/(4*Math.PI));

            let data = results_eta(a_0,n,J,g,rootprec,eps); //(a_0,n,j,g,root_prec)
            let n_s_phi = math.log(math.divide(math.multiply(math.dotMultiply(data[0],data[0]),alpha),J));
            let diff = math.abs(math.subtract(data[1],n_s_phi));
            let found = Math.min(...diff);//new spread operator
            let rho_index = diff.indexOf(found);
            let result = n_s_phi[rho_index]

            return result;
        }
        return produce_plot_floating();
    }

    function OML_calc(){  
        function OML_find_surface_potential(beta,Z){
            let m_i = Z*1.67*1e-27;
            let m_e = 9.11*1e-31;
            let rootprec = 10**(-12);
            let mu = (m_i/m_e);
            let a_0 = 1;

            let k = (((mu*beta)**0.5)*Math.exp(beta/Z))/Z;

            let W_0_H = find_W_H(a_0,k,rootprec);

            function find_eta(W_0,beta,Z){
                return W_0 - (beta/Z);
            }

            let results = find_eta(W_0_H,beta,Z);

            return results;
        }

        function OML_produce_surface_potenial_plot(){//produce data for fresnel curves
            let beta = 10**(parseFloat(document.getElementById('BetaController').value));
            let Z = parseFloat(document.getElementById('ZController').value); 
            let val = OML_find_surface_potential(beta,Z);
            return val;
        }
        return OML_produce_surface_potenial_plot();
    }

    function MOML_calc(){ 
        function MOML_find_surface_potential(beta,Z,gamma){
            let m_i = 1.67*1e-27;
            let m_e = 9.11*1e-31;
            let rootprec = 10**(-12);
            let mu = (m_i/m_e);
            let a_0 = 1;

            let c = 0.5*Math.log(2*Math.PI*(1/mu)*(1 + beta*gamma));
            let k = ((mu*beta)**0.5)*Math.exp(beta + c);

            let W_0_H = find_W_H(a_0,k,rootprec);

            function find_eta(W_0,beta,c){
                return W_0 - (beta + c);
            }

            let result = find_eta(W_0_H,beta,c);

            return result;
        }

        function MOML_produce_surface_potenial_plot(){
            let beta = 10**(parseFloat(document.getElementById('BetaController').value));
            let Z = parseFloat(document.getElementById('ZController').value); 
            let Gamma = parseFloat(document.getElementById('GammaController').value);        
            let val = MOML_find_surface_potential(beta,Z,Gamma);
            return val;
        }
        return MOML_produce_surface_potenial_plot();
    }

    function SOML_calc(){   
        function SOML_find_surface_potential(beta,Z,u){
            let m_i = 1.67*1e-27;
            let m_e = 9.11*1e-31;
            let rootprec = 10**(-12);
            let mu = (m_i/m_e);
            let a_0 = 1;

            let p_1 = (((2*Math.PI*beta)**0.5)/(4*u))*(1 + (u**2)/beta)*erf(u/((2*beta)**0.5)) + 0.5*Math.exp(-1*(u**2)/(2*beta));
            let p_2 = (((Math.PI*beta)/(2*(u**2)))**0.5)*erf(u/((2*beta)**0.5));

            let k = (((mu*beta)**0.5)/p_2)*Math.exp((beta*p_1)/p_2);

            let W_0_H = find_W_H(a_0,k,rootprec);

            function find_eta(W_0,beta,p_1,p_2){
                return W_0 - (p_1*beta/p_2);
            }
            let result = find_eta(W_0_H,beta,p_1,p_2);

            return result;
        }

        function SOML_produce_surface_potenial_plot(){//produce data for fresnel curves
            let val;
            let beta = 10**(parseFloat(document.getElementById('BetaController').value));
            let Z = parseFloat(document.getElementById('ZController').value); 
            let U =  parseFloat(document.getElementById('UController').value);
            if (U == 0){
                val = OML_find_surface_potential(beta,Z);
            }else{
                val = SOML_find_surface_potential(beta,Z,U);
            }
            return val;
        }
        return SOML_produce_surface_potenial_plot();
    }

    function SMOML_calc(){
        function SMOML_find_surface_potential(beta,Z,gamma,u){
            let m_i = 1.67*1e-27;
            let m_e = 9.11*1e-31;
            let rootprec = 10**(-12);
            let mu = (m_i/m_e);
            let a_0 = 1;

            let p_1 = (((2*Math.PI*beta)**0.5)/(4*u))*(1 + (u**2)/beta)*erf(u/((2*beta)**0.5)) + 0.5*Math.exp(-1*(u**2)/(2*beta));
            let p_2 = (((Math.PI*beta)/(2*(u**2)))**0.5)*erf(u/((2*beta)**0.5));
            let c = 0.5*Math.log(2*Math.PI*(1/mu)*(1 + beta*gamma));
            let k = (((mu*beta)**0.5)/p_2)*Math.exp(((beta*p_1)/p_2) + c);

            let W_0_H = find_W_H(a_0,k,rootprec);

            function find_eta(W_0,beta,c,p_1,p_2){
                return W_0 - (c + p_1*beta/p_2);
            }
            let result = find_eta(W_0_H,beta,c,p_1,p_2);

            return result;
        }

        function SMOML_produce_surface_potenial_plot(){
            let val;
            let beta = 10**(parseFloat(document.getElementById('BetaController').value));
            let Z = parseFloat(document.getElementById('ZController').value); 
            let Gamma = parseFloat(document.getElementById('GammaController').value);
            let U =  parseFloat(document.getElementById('UController').value);
            if (U == 0){
                val = MOML_find_surface_potential(beta,Z,Gamma);
            }else{
                val = SMOML_find_surface_potential(beta,Z,Gamma,U);
            }

            return val;
        }
        return SMOML_produce_surface_potenial_plot();
    }

    function Calculator(){
        let data;
        let selectedValue = document.getElementById("Select").value;
        switch(selectedValue) {
            case  "ABR":
                data = ABR_calc();
                break;
            case "OML":
                data = OML_calc();
                break;
            case "MOML":
                data = MOML_calc();
                break;
            case "SOML":
                data = SOML_calc();
                break;
            case "SMOML":
                data = SMOML_calc();
                break;
        }
        return data;
    }

    function update_select_sliders() {
        // NB: updates according to the active tab
        let selectedValue = document.getElementById("Select").value; // finds out which function is active
        switch(selectedValue) {
            case  "ABR":
                $('#Beta').hide();
                $('#Gamma').hide();
                $('#U').hide();
                $('#J').show();
                $('#Eps').show();
                $('#Z').show();
                break;
            case "OML":
                $('#J').hide();
                $('#Eps').hide();
                $('#Gamma').hide();
                $('#U').hide();
                $('#Beta').show();
                $('#Z').show();
                break;
            case "MOML":
                $('#J').hide();
                $('#Eps').hide();
                $('#U').hide();
                $('#Beta').show();
                $('#Z').show();
                $('#Gamma').show(); 
                break;
            case "SOML":
                $('#J').hide();
                $('#Eps').hide();
                $('#Gamma').hide();
                $('#U').show();
                $('#Beta').show();
                $('#Z').show();
                break;
            case "SMOML":
                $('#J').hide();
                $('#Eps').hide();
                $('#Beta').show();
                $('#Z').show();
                $('#Gamma').show();
                $('#U').show();
                break;
        }
    }
    function initial() {
        update_select_sliders();
        $('#Select').on("change", update_select_sliders);

        $("input[type=range]").each(function () {
            /*Allows for live update for display values*/
            $(this).on('input', function(){
                //Displays: (FLT Value) + (Corresponding Unit(if defined))
                $("#"+$(this).attr("id") + "Display").text( $(this).val() + $("#"+$(this).attr("id") + "Display").attr("data-unit"));
                //NB: Display values are restricted by their definition in the HTML to always display nice number.
            });

        });

        $('#CalcButton').on('click', function() {
            $("#Norm_surd_pot-display").html(Calculator());
        });
    }
    initial();
});