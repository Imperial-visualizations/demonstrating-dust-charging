/*jshint esversion: 6 */
$(window).on('load', function() {//main 
    let plt = {//layout of graph
        layoutW: {
            autosize: true,
            xaxis: {
                
                title: "x",
                
            },
            yaxis: {
                
                title: "W_0",
            },
            margin: {
                l: 65, r: 5, b: 40, t: 50, pad: 10
            },
        },
        layoutSurface: {
            autosize: true,
            xaxis: {
                range: [-6, 4],
                title: "Beta",
                type: 'log',
            },
            yaxis: {
                range: [0, 3],
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
            //console.log(a_n);
            a_n = a_plus;
            a_plus = calc_x_plus(a_n,z);
            //console.log(a_plus); 
        }
        W_0 = (a_n + a_plus)/2;
        return W_0; 
    }

    function find_surface_potential(beta,Z){
        let m_i = Z*1.67*1e-27;
        let m_e = 9.11*1e-31;
        let rootprec = 10**(-12);
        let mu = (m_i/m_e);
        let a_0 = 1;

        let k = (((mu*beta)**0.5)*Math.exp(beta/Z))/Z;//k becomes too larger and gets registered at infitinty beyond beta of about 340ish
        let W_0_H = find_W_H(a_0,k,rootprec);
        function find_eta(W_0,beta,Z){
            return W_0 - (beta/Z);
        }

        let result = find_eta(W_0_H,beta,Z);

        return result;
    }

    function produce_surface_potenial_plot(){//produce data for fresnel curves

        let beta = logspace(-6,4,1000);
        let z_array = numeric.linspace(1,5,5);        
        let plot_data = [];
        
        for(let j = 0;j<z_array.length;j++){
            let data_H = [];

            for(let i = 0;i<beta.length;i++){
                let val = find_surface_potential(beta[i],z_array[j]);
                data_H.push(val);
            }
            let sp_line_H = {
                x: beta,
                y: data_H,
                type: 'scatter',
                name: 'Surface potential  H Z = '+z_array[j].toString(),
            };
            plot_data.push(sp_line_H);
        }
        return plot_data;
    }
    function produce_W_plot(){
        let plot_data = [];
        let x_range = numeric.linspace(0,100,1000);
        let w_data_H = [];
        let root_prec = 10**(-8);
        let a_0 = 1;

        //console.log(upper_bound);

        for(let i = 0; i < x_range.length; i++){
            let val_H = W_0 = find_W_H(a_0,x_range[i],root_prec);
            w_data_H.push(val_H);
            //console.log(val);
        }
        let W_line_H = {
            x: x_range,
            y: w_data_H,
            type: 'scatter',
            name: 'W_0 values H',
        };
        plot_data.push(W_line_H);
        return plot_data;
    }
    function initial(){
            
        Plotly.purge("graph_W");
        Plotly.newPlot("graph_W", produce_W_plot(),plt.layoutW);

        Plotly.purge("graph_surface_potential");
        Plotly.newPlot("graph_surface_potential", produce_surface_potenial_plot(),plt.layoutSurface);

    }
    initial();
});