/*jshint esversion: 6 */
$(window).on('load', function() {//main 
    let plt = {//layout of graph
        layoutSurface: {
            autosize: true,
            xaxis: {
                //range: [0, 0.5],
                title: "Beta",
                type: 'log',
            },
            yaxis: {
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
    function find_W(x,upper_bound,root_prec){
        let root_precision = root_prec;
        
        function h(x,y){
            let val = y*Math.exp(y) - x;
            //console.log("x,y,val");
            //console.log(x,y,val);
            return val;
        }   

        function get_closer(y1,y2,x){
            let y3_gc = (y1+y2)/2;

            let s1 = Math.sign(h(x,y1));
            let s2 = Math.sign(h(x,y2));
            let s3 = Math.sign(h(x,y3_gc));

            let a;
            if (s1 == s3){
                a = y2;
            }
            else if(s2 == s3){
                a = y1;
            }
            else{
                //console.log("BROKEN");
            }
            return [a,y3_gc];
        }

        function loop(root_p, y1, y2, x){
            //console.log("y1, y2, x");
            //console.log(y1, y2, x);
            //Initial guesses x1, x2 and 'number' specifies the minimum value needed for hapiness
            let i = 0;
            while (Math.abs(y1 - y2) > root_p){
                y_res = get_closer(y1,y2,x);
                y1 = y_res[0];
                y2 = y_res[1];
                //console.log(i);
                i++;
            }
            let y3_l = (y1+y2)/2;
            return y3_l;
        }

        let y3 = loop(root_precision,0,upper_bound,x);//bisection limits set by the inital guess need to adjust this!!!!
        return y3;
    }

    function find_eta(beta,Z,mu,rootprec){
        let W_input = (mu*(beta)**0.5*Math.exp(beta/Z))/Z;
        //console.log("W_input");
        //console.log(W_input);
        let W_0 = find_W(W_input,W_input/2,rootprec);
        let eta = (beta/Z) - W_0;
        return eta;
    }

    function find_surface_potential(beta,Z){
        let m_i = Z*1.67*1e-27;
        let m_e = 9.11*1e-31;
        //let e_charge = 1.6*1e-19;
        let root_prec = 10**(-8);
        //let k_B =1.38*1e-23;
        let mu = (m_i/m_e);
        //let T_e = 1000;

        let eta_a = find_eta(beta,Z,mu,root_prec);
        //console.log("hi");
        //console.log(eta_a);
        //let phi_surface = (k_B*T_e*eta_a)/e_charge;
        //console.log(phi_surface);
        return -eta_a;
    }
    function produce_surface_potenial_plot(){//produce data for fresnel curves

        let beta = logspace(-6,4,100);
        let z_array = numeric.linspace(1,10,10);        
        let plot_data = [];
        
        for(let j = 0;j<z_array.length;j++){
            let data = [];
            for(let i = 0;i<beta.length;i++){
                data.push(find_surface_potential(beta[i],z_array[j]));
            }

            let sp_line = {
                x: beta,
                y: data,
                type: 'scatter',
                name: 'Surface potential Z = '+z_array[j].toString(),
            };
            plot_data.push(sp_line);
        }
        return plot_data;
    }
    function initial(){
            
        Plotly.purge("graph_surface_potential");
        Plotly.newPlot("graph_surface_potential", produce_surface_potenial_plot(),plt.layoutSurface);

    }
    initial();
});