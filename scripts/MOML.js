/*jshint esversion: 6 */
$(window).on('load', function() {//main 
    let plt = {//layout of graph
        layoutSurface: {
            autosize: true,
            xaxis: {
                range: [-2,1],
                title: "Beta",
                type: 'log',
            },
            yaxis: {
                range: [1, 3.5],
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

    function MOML_find_surface_potential(beta,mu,gamma){
        let rootprec = 10**(-12);
        let a_0 = 1;

        let c = 0.5*Math.log(2*Math.PI*(1/mu)*(1 + beta*gamma));
        let k = ((mu*beta)**0.5)*Math.exp(beta + c);

        let W_0_H = find_W_H(a_0,k,rootprec);

        function find_eta(W_0,beta,c){
            return W_0 - (beta + c);
        }
        
        let result = find_eta(W_0_H, beta, c);

        return result;
    }

    function MOML_produce_surface_potenial_plot(){//produce data for fresnel curves
        let gamma = numeric.linspace(0.1,3,5);
        let beta = logspace(-6,2,100);   
        let m_i = 1.67*1e-27;
        let m_e = 9.11*1e-31; 
        let mu = (m_i/m_e);     
        let plot_data = [];

        for(let j = 0; j<gamma.length ;j++){
            let data_H = [];
            for(let i = 0;i<beta.length;i++){
                let val = MOML_find_surface_potential(beta[i],mu,gamma[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: beta,
                y: data_H,
                type: 'scatter',
                name: 'MOML: Norm Surface potential H G='+ math.round(gamma[j],2).toString(),
            };
            plot_data.push(sp_line_H);
        }

        let tends_val = -0.5*math.log(2*Math.PI/mu);
        let tends_val_array = new Array (beta.length);
        tends_val_array.fill(tends_val);
        let tends = {
            x: beta,
            y: tends_val_array,
            type: 'lines',
            line: {
                dash: 'dot',
                width: 1
              },
            name: 'MOML: c',
        };
        plot_data.push(tends);
        
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

        let beta = logspace(-6,2,1000);
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

    function plot_oml_moml(){
        let oml = OML_produce_surface_potenial_plot();
        let moml = MOML_produce_surface_potenial_plot();
        return oml.concat(moml);
    }

    function initial(){

        Plotly.purge("graph_surface_potential");
        Plotly.newPlot("graph_surface_potential", plot_oml_moml(),plt.layoutSurface);

    }
    initial();
});