function logspace(a,b,c){
    let d = numeric.linspace(a,b,c)
    let e = []
    for(let i = 0;i < c;i++){
        e.push(math.pow(10,d[i]));
    }
    return e;
}
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
    function find_phi(j,g,root_prec){

        let target = j/g;
        let root_precision = root_prec;
        
        function h(x,target){
            let val = ((4*x**(3/2)*(2*x - 3)*(2*x+1)) /((2*x-1)**3)) - target;
            return val;
        }   

        function get_closer(x1,x2,target){
            let x3_gc = (x1+x2)/2;

            let s1 = Math.sign(h(x1,target));
            let s2 = Math.sign(h(x2,target));
            let s3 = Math.sign(h(x3_gc,target));

            let a;
            if (s1 == s3){
                a = x2;
            }
            else if(s2 == s3){
                a = x1;
            }
            else{
                console.log("BROKEN");
            }
            return [a,x3_gc];
        }

        function loop(root_p, x1, x2,target){
            //Initial guesses x1, x2 and 'number' specifies the minimum value needed for hapiness
            while (Math.abs(x1 - x2) > root_p){
                x_res = get_closer(x1,x2,target);
                x1 = x_res[0];
                x2 = x_res[1];
            }
            let x3_l = (x1+x2)/2;
            return x3_l;
        }

        let x3 = loop(root_precision,0.0001,0.4999,target);
        return x3;
    }

    function boundry_u(J,g,root_prec){
        let phi_calc = find_phi(J,g,root_prec);
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
    function inital_conditions(j,g,root_prec){
        let inital_u = boundry_u(j,g,root_prec);
        let inital_rho = boundry_rho(inital_u,j);
        let inital_v = boundry_v(inital_rho,inital_u,j);
        return [inital_rho, inital_u, inital_v];
    } 
    
    function f_der(rho,beta,J){
        let f = [,];
        let u = beta[0];
        let v = beta[1];
        f[0] = v;
        f[1] = rho**(-2) * J* u**(-1/2)   -   math.exp(-u)  -   2*v*rho**(-1);
        return f;
    }
    function k_1(rho,beta,h,J){
        let k_1 = math.multiply(h,f_der(rho,beta,J));
        return k_1;
    }
    function k_2(rho,beta,h,J,k_1){
        beta = math.add(beta,math.multiply(0.5,k_1));
        rho = rho+0.5*h; 
        let k_2 = math.multiply(h,f_der(rho,beta,J));
        return k_2;
    }
    function k_3(rho,beta,h,J,k_2){
        beta = math.add(beta,math.multiply(0.5,k_2));
        rho = rho+0.5*h; 
        let k_3 = math.multiply(h,f_der(rho,beta,J));
        return k_3;
    }
    function k_4(rho,beta,h,J,k_3){
        beta = math.add(beta,k_3);
        rho = rho+h; 
        let k_4 = math.multiply(h,f_der(rho,beta,J));
        return k_4;
    }
    function step(rho,u,v,h,J){
        let beta = [u,v];
        let k1 = k_1(rho,beta,h,J);
        let k2 = k_2(rho,beta,h,J,k1);
        let k3 = k_3(rho,beta,h,J,k2);
        let k4 = k_4(rho,beta,h,J,k3);

        let o_1 = math.add(k1,math.multiply(2,k2))//wont let me add all in a list, super wierd
        let o_2 = math.add(math.multiply(2,k3),k4) 
        
        beta = math.add(beta,math.multiply(1/6,math.add(o_1,o_2)));

        let new_rho = rho+h;
        let new_cond = [new_rho,beta[0],beta[1]];
        return new_cond;
    }

    function results(h,j,g,root_prec,rho_up_lim,rho_lower_lim){
        let current_step = inital_conditions(j,g,root_prec);
        let phi_data = [];
        let rho_data = [];

        if (h > 0 && current_step[0] >= rho_up_lim){
            return [phi_data,rho_data];
        }

        rho_data.push(current_step[0]);
        phi_data.push(current_step[1]);
        
        if(h > 0){
            while(current_step[0] > (rho_lower_lim + Math.abs(h)) && current_step[0] <= rho_up_lim){
                current_step = step(current_step[0],current_step[1],current_step[2],h,j);//(rho,u,v,h,J)
                rho_data.push(current_step[0]);
                phi_data.push(current_step[1]);
                
            } 
        }else{
            while(current_step[0] > (rho_lower_lim + Math.abs(h)) ){
                current_step = step(current_step[0],current_step[1],current_step[2],h,j);
                rho_data.push(current_step[0]);
                phi_data.push(current_step[1]);
            }
        }
        return [rho_data,phi_data];
    }
    function produce_floating(){
        
        let j =  parseFloat(document.getElementById('JController').value);
        let g = 10**3;
        let rootprec = 10**(-8);
        let rho_lower_lim = 0;
        let Z = 1;
        let alpha = (Z**0.5)*(1836/(4*Math.PI));
        let P_array = [];
        let h = 0.01;
        let rho_up_lim = 15;

        let data_forwards = results(h,j,g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)#
        let data_backwards = results(-h,j,g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)
        rho_data_plot_b = data_backwards[0];
        phi_data_plot_b = data_backwards[1];

        rho_data_plot_b.reverse();
        phi_data_plot_b.reverse();
        rho_data_plot_b.pop();
        phi_data_plot_b.pop();

        rho_data_plot_f = data_forwards[0];
        phi_data_plot_f = data_forwards[1];        

        let data = [rho_data_plot_b.concat(rho_data_plot_f), phi_data_plot_b.concat(phi_data_plot_f)];
        
        let n_s_phi = math.log(    math.divide(     math.multiply(      math.dotMultiply(data[0],data[0]),    alpha),      j)      );
        let diff_NaN = math.abs(math.subtract(data[1],n_s_phi));
        let diff_no_NaN  = diff_NaN.filter(function (value) {
            return !Number.isNaN(value);
        });
        let diff = math.abs(diff_no_NaN);
        let found = Math.min(...diff);//new spread operator
        let rho_index = diff.indexOf(found);
        let P_val = data[0][rho_index];
        P_array.push(P_val);
        let n_s = n_s_phi[rho_index];
        return n_s;
    }
    return produce_floating();
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
            $('#Z').hide();
            $('#Gamma').hide();
            $('#U').hide();
            $('#J').show();
            break;
        case "OML":
            $('#J').hide();
            $('#Gamma').hide();
            $('#U').hide();
            $('#Beta').show();
            $('#Z').show();
            break;
        case "MOML":
            $('#J').hide();
            $('#U').hide();
            $('#Beta').show();
            $('#Z').show();
            $('#Gamma').show(); 
            break;
        case "SOML":
            $('#J').hide();
            $('#Gamma').hide();
            $('#U').show();
            $('#Beta').show();
            $('#Z').show();
            break;
        case "SMOML":
            $('#J').hide();
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
        $("#Norm_surd_pot-display").html(Calculator().toFixed(3));
    });
}
initial();