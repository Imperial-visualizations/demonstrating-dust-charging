/*jshint esversion: 6 */
$(window).on('load', function() {//main 
    let plt = {//layout of graph
        layoutSurface: {
            autosize: true,
            xaxis: {
                range: [-6,4],
                title: "Beta",
                type: 'log',
            },
            yaxis: {
                range: [0, 4],
                title: "Normalized Surface Potential",
            },
            margin: {
                l: 65, r: 5, b: 40, t: 50, pad: 10
            },
        },
        layoutSurface_u: {
            autosize: true,
            xaxis: {
                range: [0,5],
                title: "u",
                //type: 'log',
            },
            yaxis: {
                range: [2, 4.5],
                title: "Normalized Surface Potential",
            },
            margin: {
                l: 65, r: 5, b: 40, t: 50, pad: 10
            },
        },
    };    

   function logspace(a,b,c){
        let d = numeric.linspace(a,b,c)
        let e = []
        for(let i = 0;i < c;i++){
            e.push(math.pow(10,d[i]));
        }
        return e;
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
            a_n = a_plus;
            a_plus = calc_x_plus(a_n,z);
        }
        W_0 = (a_n + a_plus)/2;
        return W_0; 
    }

    function SMOML_find_surface_potential(beta,u){
        let m_i = 1.67*1e-27;
        let m_e = 9.11*1e-31;
        let rootprec = 10**(-12);
        let mu = (m_i/m_e);
        let a_0 = 1;
        let gamma = 5/3;

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

    function SMOML_produce_surface_potenial_plot(){//produce data for fresnel curves

        let beta = logspace(-6,4,5000);
        let u_array = [1,0.5,0.25];      
        let plot_data = [];
        
        for(let j=0;j<u_array.length;j++){
            let data_H = [];
            for(let i = 0;i<beta.length;i++){
                let val = SMOML_find_surface_potential(beta[i],u_array[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: beta,
                y: data_H,
                type: 'scatter',
                name: 'SMOML: Norm Surface potential H, u = ' + u_array[j].toString(),
            };
            plot_data.push(sp_line_H);
        }
        return plot_data;
    }

    function SOML_find_surface_potential(beta,u){
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

        let beta = logspace(-6,4,5000);
        let u_array = [1,0.5,0.25];      
        let plot_data = [];
        
        for(let j=0;j<u_array.length;j++){
            let data_H = [];
            for(let i = 0;i<beta.length;i++){
                let val = SOML_find_surface_potential(beta[i],u_array[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: beta,
                y: data_H,
                type: 'scatter',
                name: 'SOML: Norm Surface potential H, u = ' + u_array[j].toString(),
            };
            plot_data.push(sp_line_H);
        }
        return plot_data;
    }

    function MOML_find_surface_potential(beta){
        let m_i = 1.67*1e-27;
        let m_e = 9.11*1e-31;
        let rootprec = 10**(-12);
        let gamma = 5/3;
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

    function MOML_produce_surface_potenial_plot(){//produce data for fresnel curves

        let beta = logspace(-6,4,5000);        
        let plot_data = [];
        let data_H = [];

        for(let i = 0;i<beta.length;i++){
            let val = MOML_find_surface_potential(beta[i]);
            data_H.push(val);
        }

        let sp_line_H = {
            x: beta,
            y: data_H,
            type: 'scatter',
            name: 'MOML: Norm Surface potential H',
        };

        plot_data.push(sp_line_H);
        
        return plot_data;
    }

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

        let beta = logspace(-6,4,5000);
        let z = 1;      
        let plot_data = [];
        let data_H = [];

        for(let i = 0;i<beta.length;i++){
            let val = OML_find_surface_potential(beta[i],z);
            data_H.push(val);
        }

        let sp_line_H = {
            x: beta,
            y: data_H,
            type: 'scatter',
            name: 'OML: Norm Surface potential H',
        };

        plot_data.push(sp_line_H);
        return plot_data;
    }

    function SMOML_produce_surface_potenial_plot_u(){//produce data for fresnel curves

        let beta = [1,10];
        let u_array = numeric.linspace(0,5,1000);      
        let plot_data = [];
        
        for(let i = 0;i<beta.length;i++){
            let data_H = [];
            for(let j=0;j<u_array.length;j++){
                let val = SMOML_find_surface_potential(beta[i],u_array[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: u_array,
                y: data_H,
                type: 'scatter',
                name: 'SMOML: Norm Surface potential H, beta = ' + beta[i].toString(),
            };
            plot_data.push(sp_line_H);
        }
        return plot_data;
    }
    function SOML_produce_surface_potenial_plot_u(){//produce data for fresnel curves

        let beta = [1,10];
        let u_array = numeric.linspace(0,5,1000);       
        let plot_data = [];
        
        for(let i = 0;i<beta.length;i++){
            let data_H = [];
            for(let j=0;j<u_array.length;j++){
                let val = SOML_find_surface_potential(beta[i],u_array[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: u_array,
                y: data_H,
                type: 'scatter',
                name: 'SOML: Norm Surface potential H, beta = ' + beta[i].toString(),
            };
            plot_data.push(sp_line_H);
        }
        return plot_data;
    }


    function plot_oml_moml_soml_smoml(){
        let oml = OML_produce_surface_potenial_plot();
        let moml = MOML_produce_surface_potenial_plot();
        let soml = SOML_produce_surface_potenial_plot();
        let smoml = SMOML_produce_surface_potenial_plot();
        let data_surface = oml.concat(moml,soml,smoml);
        return data_surface;
    }

    function plot_soml_smoml_u(){
        let soml = SOML_produce_surface_potenial_plot_u();
        let smoml = SMOML_produce_surface_potenial_plot_u();
        let data_surface = soml.concat(smoml);
        return data_surface;
    }

    function initial(){

        Plotly.purge("graph_surface_potential");
        Plotly.newPlot("graph_surface_potential", plot_oml_moml_soml_smoml(),plt.layoutSurface);

        Plotly.purge("graph_surface_potential_u");
        Plotly.newPlot("graph_surface_potential_u", plot_soml_smoml_u(),plt.layoutSurface_u);

    }
    initial();
});