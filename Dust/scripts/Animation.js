
/*jshint esversion:7*/
//set global variables
//allpoints for storing charges, maxpoints to limit total n of allpoints, newchargex/y for position of new charge on top

let width = $('#sketch-holder').width(), height = $('#sketch-holder').height(), allpoints = [], newchargex = width*0.96, newchargey = height*0.27, spam_num = 20, source_num = 20, Bx = 0, By = 0,mouse_listx = [],mouse_listy = [],trace_freq = 3,trace_freq_crit = 2,source_x = 100;
const R = 4, rect_height = height*0.3, grain_x = width/2, grain_y = (height + rect_height)/2,grain_z = 0, grain_R = 50,trace_len = 100,mouse_len = 5, mag_repx = 50, mag_repy = rect_height + 50;

const m_i = 1.67*1e-27;
const m_e = 9.11*1e-31;
const mu = (m_i/m_e);
const epsilon_0 = 8.85*1e-12;
const k_b = 1.38*1e-23;

const B_field_adjust = 1e-7;
const color_shift = 1e10;

const vel_factor = 1/2;
const e_charge = -1.6*1e-19;
const i_charge = 1.6*1.e-19;
const n_0 = 1;

let T_e = 1e2;
const lambda_de = Math.sqrt((epsilon_0*k_b*T_e)/(n_0*(e_charge**2)));
let T_i = T_e*10**(parseFloat(document.getElementById('BetaController').value));
let lambda_di = Math.sqrt((epsilon_0*k_b*T_i)/(n_0*(e_charge**2)));
let lambda_d = Math.sqrt(1/(1/(lambda_de**2) + 1/(lambda_di**2)));

const root =  10**(-10);
const eps = 10**(-10); 
let dt = 10**parseFloat(document.getElementById('PlaySpeed').value);

let Bz;
let Q_grain;
let Grain_colour;
let phi_grain;
let q_charx_list_1;
let q_chary_list_1;
let q_charx_list_2;
let q_chary_list_2;

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
            return [beta_5,beta_4];
        }

        function step(rho,u,v,h,J,eps){
            let delta_0 = math.abs(math.multiply(eps,[u,v]));
            let S = 0.95;
            let h_new;
            let rk = RK_fomulae(rho,u,v,h,J);
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
        let a_0 = 0.1;
        let n = 5000;        
        //let eps = 10**parseFloat(document.getElementById('EpsController').value); 
        let J = 10**parseFloat(document.getElementById('JController').value); 
        let Z = parseFloat(document.getElementById('ZController').value); 
        //let root =  10**(parseFloat(document.getElementById('RootController').value));
        
        let alpha = (Z**0.5)*(1836/(4*Math.PI));

        let data = results_eta(a_0,n,J,g,root,eps); //(a_0,n,j,g,root_prec)
        let n_s_phi = math.log(math.divide(math.multiply(math.dotMultiply(data[0],data[0]),alpha),J));
        let diff = math.abs(math.subtract(data[1],n_s_phi));
        let found = Math.min(...diff);//new spread operator
        let rho_index = diff.indexOf(found);
        let result = n_s_phi[rho_index]

        return result;
    }
    return produce_plot_floating();
}

function OML_find_surface_potential(beta,Z,mu,rootprec){
    let a_0 = 1;

    let k = (((mu*beta)**0.5)/Z)*Math.exp(beta/Z);

    let W_0_H = find_W_H(a_0,k,rootprec);

    let results = W_0_H - (beta/Z);

    return results;
}

function OML_produce_surface_potenial_plot(){
    
    let beta = 10**(parseFloat(document.getElementById('BetaController').value));
    let Z = 1;
    
    let val = OML_find_surface_potential(beta,Z,mu,root);
    return val;
}

function MOML_find_surface_potential(beta,Z,gamma,mu,rootprec){
    let a_0 = 1;

    let c = 0.5*Math.log(2*Math.PI*(1/mu)*(1 + beta*gamma));
    let k = (((mu*beta)**0.5)/Z)*Math.exp(beta/Z + c);

    let W_0_H = find_W_H(a_0,k,rootprec);
    let result =  W_0_H - (beta/Z + c);

    return result;
}

function MOML_produce_surface_potenial_plot(){
    let beta = 10**(parseFloat(document.getElementById('BetaController').value));
    let Z = 1;
    let Gamma = parseFloat(document.getElementById('GammaController').value);   

    let val = MOML_find_surface_potential(beta,Z,Gamma,mu,root);
    return val;
}


function SOML_find_surface_potential(beta,Z,u,mu,rootprec){
    let a_0 = 1;

    let p_1 = (((2*Math.PI*beta)**0.5)/(4*u))*(1 + (u**2)/beta)*erf(u/((2*beta)**0.5)) + 0.5*Math.exp(-1*(u**2)/(2*beta));
    let p_2 = (((Math.PI*beta)/(2*(u**2)))**0.5)*erf(u/((2*beta)**0.5));

    let k = (((mu*beta)**0.5)/(Z*p_2))*Math.exp((beta*p_1)/(Z*p_2));

    let W_0_H = find_W_H(a_0,k,rootprec);

    let result = W_0_H - ((p_1*beta)/(Z*p_2));

    return result;
}

function SOML_produce_surface_potenial_plot(){//produce data for fresnel curves
    let val;

    let beta = 10**(parseFloat(document.getElementById('BetaController').value));
    let Z = 1;
    let U =  parseFloat(document.getElementById('UController').value);
    //let root =  10**(parseFloat(document.getElementById('RootController').value));

    if (U == 0){
        val = OML_find_surface_potential(beta,Z,mu,root);
    }else{
        val = SOML_find_surface_potential(beta,Z,U,mu,root);
    }
    return val;
}

function SMOML_find_surface_potential(beta,Z,gamma,u,mu,rootprec){
    let a_0 = 1;

    let p_1 = (((2*Math.PI*beta)**0.5)/(4*u))*(1 + (u**2)/beta)*erf(u/((2*beta)**0.5)) + 0.5*Math.exp(-1*(u**2)/(2*beta));
    let p_2 = (((Math.PI*beta)/(2*(u**2)))**0.5)*erf(u/((2*beta)**0.5));
    let c = 0.5*Math.log(2*Math.PI*(1/mu)*(1 + beta*gamma));
    let k = (((mu*beta)**0.5)/(Z*p_2))*Math.exp(((beta*p_1)/(Z*p_2)) + c);

    let W_0_H = find_W_H(a_0,k,rootprec);

    let result = W_0_H - (c + ((p_1*beta)/(Z*p_2)));;

    return result;
}

function SMOML_produce_surface_potenial_plot(){
    let val;
    let beta = 10**(parseFloat(document.getElementById('BetaController').value));
    let Z = 1;
    let Gamma = parseFloat(document.getElementById('GammaController').value);
    let U =  parseFloat(document.getElementById('UController').value);

    if (U == 0){
        val = MOML_find_surface_potential(beta,Z,Gamma,mu,root);
    }else{
        val = SMOML_find_surface_potential(beta,Z,Gamma,U,mu,root);
    }

    return val;
}

function Calculator(){
    let data;
    let selectedValue_model = document.getElementById("Select_model").value;
    //console.log(selectedValue_model);
    switch(selectedValue_model) {
        case  "ABR":
            data = ABR_calc();
            break;
        case "OML":
            data = OML_produce_surface_potenial_plot();
            break;
        case "MOML":
            data = MOML_produce_surface_potenial_plot();
            break;
        case "SOML":
            data = SOML_produce_surface_potenial_plot();
            break;
        case "SMOML":
            data = SMOML_produce_surface_potenial_plot();
            break;
    }
    return data;
}

function update_select_sliders() {
    // NB: updates according to the active tab
    $("#grain_charge").html();
    let selectedValue = document.getElementById("Select_model").value; // finds out which function is active
    switch(selectedValue) {
        case  "ABR":
            $('#Beta').hide();
            $('#Gamma').hide();
            $('#U').hide();
            $('#J').show();
            $('#Z').show();
            
            break;
        case "OML":
            $('#J').hide();
            $('#Gamma').hide();
            $('#U').hide();
            $('#Z').hide();
            $('#Beta').show();
            
            break;
        case "MOML":
            $('#J').hide();
            $('#U').hide();
            $('#Z').hide();
            $('#Beta').show();
            $('#Gamma').show(); 
            
            break;
        case "SOML":
            $('#J').hide();
            $('#Gamma').hide();
            $('#Z').hide();
            $('#U').show();
            $('#Beta').show();
            
            break;
        case "SMOML":
            $('#J').hide();
            $('#Z').hide();
            $('#Beta').show();
            $('#Gamma').show();
            $('#U').show();
            
            break;
    }
}

function grain_colour(){
    if (Q_grain == 0){
        this.color = "#00FF00";
    } else if (Q_grain > 0){
        let tune1 = Math.round(180 - 120*(1-Math.exp(-Math.abs(Q_grain))));
        let tune2 = Math.round(90*(Math.exp(-Math.abs(Q_grain))));
        Grain_colour = "rgb(255," + tune1.toString() + "," + tune2.toString() + ")";
    } else {
        let tune1 = Math.round(100 - 100*Math.exp(Q_grain*color_shift));
        let tune2 = Math.round(200 - 40*Math.exp(Q_grain*color_shift));
        Grain_colour = "rgb(" + tune1.toString() + "," + tune2.toString() + ",255)";
        //console.log(Grain_colour);
    }
}

function calc_charge_grain(){
    let n_float = Calculator();
    phi_grain = (k_b*T_e*n_float)/(e_charge);
    Q_grain = 4*Math.PI*epsilon_0*grain_R*phi_grain;
    $("#grain_charge").html(Q_grain.toExponential(3));
    grain_colour();
}

class charge {
    constructor(q, x, y, z, vx, vy , vz, m){

        this.positionsx = [];
        this.positionsy = [];
        this.trace_index = 0;

        this.q = q;
        this.x = x;
        this.y = y;
        this.z = z;
        this.radial = Math.sqrt(((grain_x - this.x)**2 + (grain_y - this.y)**2 + (grain_z - this.z)**2));
        this.vx = vx;
        this.vy = vy;
        this.vz = vz;
        this.m = m;
        this.r = R;
        this.clicked = false;

        if (q > 0){
            let tune1 = Math.round(180 - 120*(1-Math.exp(-Math.abs(q))));
            let tune2 = Math.round(90*(Math.exp(-Math.abs(q))));
            this.color = "rgb(255," + tune1.toString() + "," + tune2.toString() + ")";
        } else if (q < 0){
            let tune1 = Math.round(120*(Math.exp(-Math.abs(q))));
            let tune2 = Math.round((180 - 120*(1-Math.exp(-Math.abs(q)))));
            this.color = "rgb(" + tune1.toString() + "," + tune2.toString() + ",255)";
        } else {
            this.color = "#00FF00";
        }
    }

    //Cursor interactivity
    pressed(){
        if (dist(mouseX, mouseY, this.x, this.y) < this.r){
            this.clicked = true;
        }
    }

    dragposition(){
            this.x = mouseX;
            this.y = mouseY;
    }

    coulomb_force(x,y,z){
        let r = Math.sqrt(((grain_x - x)**2 + (grain_y - y)**2 + (grain_z - z)**2));
        let Fx = (1/(4*Math.PI*epsilon_0))*(Q_grain*this.q)*(x - grain_x) / (Math.pow(r,3)); 
        let Fy = (1/(4*Math.PI*epsilon_0))*(Q_grain*this.q)*(y - grain_y) / (Math.pow(r,3)); 
        let Fz = (1/(4*Math.PI*epsilon_0))*(Q_grain*this.q)*(z - grain_z) / (Math.pow(r,3)); 
        return[Fx,Fy,Fz];
    }
    deybe_huckle(x,y,z){
        let r = Math.sqrt(((grain_x - x)**2 + (grain_y - y)**2 + (grain_z - z)**2));
        let Fx = this.q*phi_grain*grain_R*Math.exp((grain_R/lambda_d)-(r/lambda_d))*((1/(r**2)) + 1/(lambda_d*r))*((x - grain_x)/r); 
        let Fy = this.q*phi_grain*grain_R*Math.exp((grain_R/lambda_d)-(r/lambda_d))*((1/(r**2)) + 1/(lambda_d*r))*((y - grain_y)/r); 
        let Fz = this.q*phi_grain*grain_R*Math.exp((grain_R/lambda_d)-(r/lambda_d))*((1/(r**2)) + 1/(lambda_d*r))*((z - grain_z)/r);  
        return[Fx,Fy,Fz];
    }

    Bfield_cross(vx,vy,vz){
        let v = [vx,vy,vz];
        let B = [Bx,By,Bz];
        let F = math.multiply(this.q,math.cross(v,B));
        return F;
    }

    f_der(W_vec){
        let x = W_vec[0];
        let y = W_vec[1];
        let z = W_vec[2];
        let vx = W_vec[3];
        let vy = W_vec[4];
        let vz = W_vec[5];

        let f = [vx,vy,vz];
        let selectedValue_field = document.getElementById("Select_field").value;
        if(selectedValue_field === "Coulomb"){
            if (document.getElementById('BOption').checked == true && Bz != 0) {
                f = f.concat(math.multiply(1/this.m,math.add(this.coulomb_force(x,y,z),this.Bfield_cross(vx,vy,vz))));
            }
            else{
                f = f.concat(math.multiply(1/this.m,this.coulomb_force(x,y,z)));
            }
        }else{
            if (document.getElementById('BOption').checked == true && Bz != 0) {
                f = f.concat(math.multiply(1/this.m,math.add(this.deybe_huckle(x,y,z),this.Bfield_cross(vx,vy,vz))));
            }
            else{
                f = f.concat(math.multiply(1/this.m,this.deybe_huckle(x,y,z)));
            }
        }
        return f;
    }
    k_1(W_vec){
        let res = this.f_der(W_vec);
        let k_1 = math.multiply(dt,res);
        return k_1;
    }
    k_2(W_vec,k_1){
        W_vec = math.add(W_vec,math.divide(k_1,2));
        let res = this.f_der(W_vec);
        let k_2 = math.multiply(dt,res);
        return k_2;
    }
    k_3(W_vec,k_2){
        W_vec = math.add(W_vec,math.divide(k_2,2));
        let res = this.f_der(W_vec);
        let k_3 = math.multiply(dt,res);
        return k_3;
    }
    k_4(W_vec,k_3){
        W_vec = math.add(W_vec,k_3);
        let res = this.f_der(W_vec);
        let k_4 = math.multiply(dt,res);
        return k_4;
    }

    step(x,y,z,vx,vy,vz){
        let W_vec = [x,y,z,vx,vy,vz];
        let k1 = this.k_1(W_vec);
        let k2 = this.k_2(W_vec,k1);
        let k3 = this.k_3(W_vec,k2);
        let k4 = this.k_4(W_vec,k3);
        W_vec = math.add(W_vec,math.multiply(1/6,math.add(k1,math.multiply(2,k2),math.multiply(2,k3),k4)));
        return W_vec;
    }
    
    motion(){
        let new_W_vec = this.step(this.x,this.y,this.z,this.vx,this.vy,this.vz);
        this.x = new_W_vec[0];
        this.y = new_W_vec[1];
        this.z = new_W_vec[2];
        this.vx = new_W_vec[3];
        this.vy = new_W_vec[4];
        this.vz = new_W_vec[5];
        this.radial = Math.sqrt(((grain_x - this.x)**2 + (grain_y - this.y)**2 + (grain_z - this.z)**2));

        if (document.getElementById('TraceOption').checked == true) {
            this.trace_index += 1;
            if (this.positionsx.length === trace_len){
                this.positionsx.shift();
                this.positionsy.shift();            
            }
            if (this.trace_index%trace_freq === 0){
                this.positionsx.push(this.x);
                this.positionsy.push(this.y);
            }
        }
    }
}


//Selects the charge that user wants
class charge_selector{
    constructor(q, x, y, z){

        this.q = q;
        this.x = x;
        this.y = y;
        this.z = z;
        this.r = R;
        this.clicked = false;

        //Colour of charge in relation to magnitude and polarity
        if (q == 0){
            this.color = "#00FF00";
        } else if (q > 0){
            let tune1 = Math.round(180 - 120*(1-Math.exp(-Math.abs(q))));
            let tune2 = Math.round(90*(Math.exp(-Math.abs(q))));
            this.color = "rgb(255," + tune1.toString() + "," + tune2.toString() + ")";
        } else {
            let tune1 = Math.round(120*(Math.exp(-Math.abs(q))));
            let tune2 = Math.round((180 - 120*(1-Math.exp(-Math.abs(q)))));
            this.color = "rgb(" + tune1.toString() + "," + tune2.toString() + ",255)";
        }
    }
    
    //God's hand to pull charge out of thin air 
    pressed(){
        if (dist(mouseX, mouseY, this.x, this.y) < this.r){
            let q_mass;
            if (parseFloat(document.getElementById('magnit').value) < 0){
                q_mass = m_e;
            }
            else{
                q_mass = m_i;
            }
            let q = new charge(this.q, this.x, this.y, this.z,0,0,0,q_mass);
            allpoints.push(q);
        }
    }
}

//function that 'move' a charge when it is clicked
function mousePressed(){
    sel.pressed();
    for (let i = 0; i < allpoints.length; i++) {
        allpoints[i].pressed();
    }
}

function mouseReleased() {
    for (let i = 0; i < allpoints.length; i++) {
        if (allpoints[i].y < rect_height || allpoints[i].y > height|| allpoints[i].x > width || allpoints[i].x < 0 ){
            allpoints.splice(i,1);
        } else {
            allpoints[i].clicked = false;
        }
    }
}

function spam_produce(){
    let vx_charge,vy_charge;   
    let x_charge = Math.random()*width;
    let y_charge = Math.random()*height;
    let therm_speed_e = Math.sqrt((k_b*T_e)/m_e);
    let therm_speed_i = Math.sqrt((k_b*T_i)/m_i);
    for(let i = 0 ; i < spam_num;i++){//q, x, y, vx, vy , m
        let charge_q =  math.random(i_charge*parseFloat(document.getElementById('magnit').min),i_charge*parseFloat(document.getElementById('magnit').max));
        let q_mass;
        if (document.getElementById("Select_model").value === "ABR"){
            vx_charge = 0;
            vy_charge = 0;
            if (charge_q < 0){
                q_mass = m_e;
            }
            else{
                q_mass = m_i;
            }
        }else{
            if(charge_q < 0){
                vx_charge = math.random(-therm_speed_e,therm_speed_e); 
                vy_charge = math.random(-therm_speed_e,therm_speed_e); 
                q_mass = m_e;
             }
             else{
                vx_charge = math.random(-therm_speed_i,therm_speed_i); 
                vy_charge = math.random(-therm_speed_i,therm_speed_i); 
                 q_mass = m_i; 
             }
        }
        let q = new charge(charge_q ,x_charge,y_charge, 0 ,vx_charge,vy_charge, 0, q_mass);
        console.log(charge_q ,x_charge,y_charge, 0 ,vx_charge,vy_charge, 0, q_mass);
        allpoints.push(q);        
    }
}

function source_produce(){
    let vx_charge;  
    let x_charge = source_x;
    let y_charge = numeric.linspace(rect_height,height,source_num);
    let z_charge = 0;
    let vz_charge = 0;
    let vy_charge = 0;
    let charge_q = i_charge*parseFloat(document.getElementById('magnit').value);

    let q_mass;
    if (charge_q < 0){
        q_mass = m_e;
    }
    else{
        q_mass = m_i;
    }
    for(let i = 0 ; i < y_charge.length;i++){//q, x, y, vx, vy , m
        if (document.getElementById("Select_model").value === "ABR"){
            vx_charge = 0;
        }else{
            if(charge_q < 0){
               vx_charge = Math.sqrt((k_b*T_e)/m_e); 
            }
            else{
                vx_charge = Math.sqrt((k_b*T_i)/m_i); 
            }
        }
        let q = new charge(charge_q ,x_charge,y_charge[i],z_charge,vx_charge,vy_charge, vz_charge,q_mass);
        allpoints.push(q);        
    }
}

function bCrit_array(){
    let vy_charge = 0,vz_charge = 0,z_charge = 0;
    let v = Math.sqrt((k_b*T_i)/m_i);
    
    let charge_q = i_charge*parseFloat(document.getElementById('magnit').value);
    let b_crit = (grain_R/v)*Math.sqrt((v**2) - ((2*charge_q*phi_grain)/m_i));
    
    let q_1 = new charge(charge_q,source_x,grain_y - b_crit ,z_charge,v,vy_charge,vz_charge,m_i);
    let q_2 = new charge(charge_q,source_x,grain_y + b_crit ,z_charge,v,vy_charge,vz_charge,m_i);

    let r_list_1 = [q_1.radial];
    q_charx_list_1 = [q_1.x];
    q_chary_list_1 = [q_1.y];

    q_charx_list_2 = [q_2.x];
    q_chary_list_2 = [q_2.y];

    let b_trace_index = 0;

    while((Math.sqrt(((grain_x - q_1.x) ** 2 + (grain_y - q_1.y) ** 2 + (grain_z - q_1.z) ** 2)) >= 1.01*(q_1.r + grain_R))) { 
        q_1.motion();
        q_2.motion();
        b_trace_index += 1;
        r_list_1.push(q_1.radial- (q_1.r + grain_R));
        if (b_trace_index%trace_freq_crit === 0){
            q_charx_list_1.push(q_1.x);
            q_chary_list_1.push(q_1.y);
            q_charx_list_2.push(q_2.x);
            q_chary_list_2.push(q_2.y);
        }
        if(r_list_1[r_list_1.length-1] - r_list_1[r_list_1.length-2] > 0){
            break;
        }
    }
}

function clear_pos(){
    for (let i = 0; i < allpoints.length; i++) {
        allpoints[i].positionsx = [];
        allpoints[i].positionsy = [];
    }
}

function mouseflick_vel(){
    let vx = 0,vy = 0;
    if(document.getElementById("Select_model").value === "ABR"){
        vx = 0;
        vy = 0;
    }else{
        for(let i = 0; i < (mouse_listx.length -1);i++){
            vx += mouse_listx[i+1] - mouse_listx[i];
            vy += mouse_listy[i+1] - mouse_listy[i];
        }
        vx = ((vx/(dt))/mouse_listx.length)*vel_factor;
        vy = ((vy/(dt))/mouse_listx.length)*vel_factor;        
    }
    return [vx,vy];
}

class currentLoop {
    constructor(I, x, y) {
        this.In = I;
        this.x = x;
        this.y = y;
        this.a = 50;
    }

    drawLoop() {
        if (this.In > 0) {

            push();
            noFill();
            stroke(212, 36, 36);
            strokeWeight(2);
            //ellipse(this.x, this.y, this.a+40, this.a+40)
            ellipse(this.x, this.y, this.a, this.a);
            pop();
            push(); // Start a new drawing state
            noStroke();
            fill(212, 36, 36);
            strokeWeight(2);
            let b = Math.round(this.a / 3);
            ellipse(this.x, this.y, b, b);
            pop();

        } else if(this.In < 0) {
            push();
            strokeWeight(2);
            noFill();
            stroke(44, 61, 171);
            //ellipse(this.x, this.y, this.a+40, this.a+40)
            ellipse(this.x, this.y, this.a, this.a);
            let c = this.a * Math.cos(Math.PI / 4) / 2;
            line(this.x + c, this.y + c, this.x - c, this.y - c);
            line(this.x - c, this.y + c, this.x + c, this.y - c);

            strokeWeight(2);
            pop();
        } else {
            push();
            noFill();
            stroke(114, 212, 119);
            strokeWeight(2);
            ellipse(this.x, this.y, this.a, this.a);
            pop();
        }
    }

    getPos() {
        return [this.x, this.y]
    }
}

//draw canvas in which everything p5.js happens
function setup() {
    update_select_sliders();
    calc_charge_grain();
    
    
    let canvas = createCanvas(width,height);
    canvas.parent('sketch-holder');
    frameRate(60);

    $("#BetaController").on("change",function(){
        T_i = T_e*10**(parseFloat(document.getElementById('BetaController').value));
        lambda_di = Math.sqrt((epsilon_0*k_b*T_i)/(n_0*(e_charge**2)));
        lambda_d = Math.sqrt(1/(1/(lambda_de**2) + 1/(lambda_di**2)));
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }
        });
    $("#Select_model").on("change",update_select_sliders);
    $("#calcbutton").on("click",function(){
        calc_charge_grain();
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }

    });
    $("#magnit").on("change",function(){
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }});
    $("#mass").on("change",function(){
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }});
    $("#vxspeed").on("change",function(){
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }});
    $("#CritOption").on("change",function(){
        if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR"){
            bCrit_array();
        }});
    $("#PlaySpeed").on("change",function(){
        dt = 10**parseFloat(document.getElementById('PlaySpeed').value);
        });
    $("#spambutton").on("click",spam_produce);
    $("#sourcebutton").on("click",source_produce);
    $("#resetbutton").on("click",function(){
        allpoints = [];
    });
    $("#TraceOption").on("click", clear_pos);
}

//main function that repeats as soon as the last line is called
function draw() {
    clear();

    if(document.getElementById('BOption').checked == true){
        Bz = parseFloat(document.getElementById('Bfield').value)*B_field_adjust;
        let mag_rep = new currentLoop(Bz,mag_repx,mag_repy);
        mag_rep.drawLoop();        
    }

    //DRAW GRAIN
    fill(color(Grain_colour));
    ellipse(grain_x, grain_y, 2*grain_R);

    if(document.getElementById('CritOption').checked == true && document.getElementById("Select_model").value != "ABR" && document.getElementById('BOption').checked == false){
        for(let j = 0;j < q_charx_list_1.length;j++){
            fill(7, 8, 8);
            ellipse(q_charx_list_1[j], q_chary_list_1[j], 2);
            ellipse(q_charx_list_2[j], q_chary_list_2[j], 2);
        }
    }

    //Brings in user input and turn it into a charge
    sel = new charge_selector(i_charge*parseFloat(document.getElementById('magnit').value), newchargex, newchargey,0);

    for (let i = 0; i < allpoints.length; i++) {
        if (allpoints[i].clicked == true){
            allpoints[i].dragposition();
            if (mouse_listx.length === mouse_len){
                mouse_listx.shift();
                mouse_listy.shift();
            }
            mouse_listx.push(mouseX);
            mouse_listy.push(mouseY);
            let mouse_speed = mouseflick_vel();
            allpoints[i].vx = mouse_speed[0];
            allpoints[i].vy = mouse_speed[1];
        }
        else{
            allpoints[i].motion();
            if (Math.sqrt(((grain_x - allpoints[i].x) ** 2 + (grain_y - allpoints[i].y) ** 2 + (grain_z - allpoints[i].z) ** 2)) <= (allpoints[i].r + grain_R) || allpoints[i].x >= width || allpoints[i].x <= 0 || allpoints[i].y <= rect_height || allpoints[i].y >= height){
                allpoints.splice(i,1);
                continue;
            }
            if (document.getElementById('TraceOption').checked == true) {
                for(let j = 0;j < allpoints[i].positionsx.length;j++){
                    fill(20, 200, 100);
                    ellipse(allpoints[i].positionsx[j], allpoints[i].positionsy[j], 2);
                }
            }
        }
        noStroke(1);
        fill(color(allpoints[i].color));
        ellipse(allpoints[i].x, allpoints[i].y, R*2);
    }

    //draw the top blue rectangle box that contains the text, slider and new charge
    noStroke();
    fill(213, 219, 222);
    rect(0, 0, width, rect_height);

    stroke(72, 99, 95);
    line(0, rect_height, width, rect_height);

    //draw the charges that are already inside the canvas

    noStroke();
    fill(color(sel.color));
    ellipse(sel.x, sel.y, R*2);
    
}

