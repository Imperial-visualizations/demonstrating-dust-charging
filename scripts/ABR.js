/*jshint esversion: 6 */
$(window).on('load', function() {//main  

    const dom = {
        gSlider: $("input#gamma"),//gamma slider
        jSlider: $("input#jnorm"),//j slider
    };

    let plt = {//layout of graph
            layoutJ: {
                autosize: true,
                xaxis: {
                    range: [0, 0.5],
                    title: "Phi_b",
                },
                yaxis: {
                    type: 'log',
                    title: "J/gamma",
                    range: [-10, 10],
                },
                margin: {
                    l: 65, r: 5, b: 40, t: 50, pad: 10
                },
            },
            layoutBray: {
                autosize: true,
                xaxis: {
                    range: [0, 15],
                    title: "Rho",
                },
                yaxis: {
                    title: "Phi",
                    range: [0, 10],
                },
                margin: {
                    l: 40, r: 5, b: 40, t: 50, pad: 10
                },
            },
            layoutKA: {
                autosize: true,
                xaxis: {
                    range: [0, 8],
                    title: "Rho",
                },
                yaxis: {
                    title: "Phi",
                    range: [0, 5],
                },
                margin: {
                    l: 40, r: 5, b: 40, t: 50, pad: 10
                },
            },
            layoutFloating: {
                autosize: true,
                xaxis: {
                    title: "P",
                    type: 'log',
                },
                yaxis: {
                    title: "Phi Surface",
                },
                margin: {
                    l: 40, r: 5, b: 40, t: 50, pad: 10
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

        //console.log("beta uv");
        //console.log(beta);
        beta = math.add(beta,math.multiply(1/6,math.add(k1,math.multiply(2,k2),math.multiply(2,k3),k4)));
        //console.log(beta);
        let new_rho = rho+h;
        let new_cond = [new_rho,beta[0],beta[1]];
        return new_cond;
    }

    function results(h,j,g,root_prec,rho_up_lim,rho_lower_lim){
        //console.log(j);
        let current_step = inital_conditions(j,g,root_prec);
        //console.log("current_step - [inital_rho,inital_u,inital_v]");
        //console.log(current_step);
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

    function produce_plot(){//produce data for fresnel curves
        $("#gamma-display").html($("input#gamma").val());//update display
        let plot_data = [];
        //[0.1,5,30,55,80] bray
        //[0.1,0.3,1,3,10] KA

        let g = 10**(parseFloat($("input#gamma").val()));

        let j_range = [0.1,0.3,1,3,10];
        let rootprec = 10**(-8);
        let h = 0.1;
        let rho_up_lim = 15
        let rho_lower_lim = 0
        for (let i = 0; i< j_range.length; i++){
            //console.log(j_range[i]);
            let data_forwards = results(h,j_range[i],g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)
            let data_backwards = results(-h,j_range[i],g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)
            //console.log(data_forwards);
            //console.log(data_backwards);

            let init_val_rho = data_backwards[0][0];
            let init_val_phi = data_backwards[1][0];
            
            //console.log("rho,phi");
            //console.log(init_val_rho,init_val_phi);

            rho_data_plot_b = data_backwards[0];
            phi_data_plot_b = data_backwards[1];

            rho_data_plot_b.reverse();
            phi_data_plot_b.reverse();
            rho_data_plot_b.pop();
            phi_data_plot_b.pop();

            rho_data_plot_f = data_forwards[0];
            phi_data_plot_f = data_forwards[1];
            //console.log(rho_data_plot_b);
            //console.log(rho_data_plot_f);
            
            let data = [rho_data_plot_b.concat(rho_data_plot_f), phi_data_plot_b.concat(phi_data_plot_f)];
            
            //console.log(data);

            let v_line = {
                x: data[0],
                y: data[1],
                type: 'scatter',
                name: 'Phi curve J = ' + j_range[i].toString(),
            };
            let marker_init = {//marker follows inital value
                x: [init_val_rho],
                y: [init_val_phi],
                opacity: 0.5,
                showlegend: false,
                type: "scatter",
                mode:"markers",
                name: 'Boundry condition'+ j_range[i].toString(),
                marker: {color: "#002147", size: 6}
            };
            //console.log(init_val[2],init_val[1]);

            plot_data.push(v_line);
            plot_data.push(marker_init);
        };
        //console.log(plot_data);
        return plot_data;
    };

    function produce_J_plot(){//produce data for fresnel curves

        $("#gamma-display").html($("input#gamma").val());//update display

        let g = 10**(parseFloat($("input#gamma").val()));
        //console.log(g);
        let plot_data = [];
        let data_eq = [];
        let data_root = [];
        let phi_range = numeric.linspace(0,0.5,1000);
        let root_prec = 10**(-8);

        function equation_j_g(x){
            j_g =  (4*x**(3/2)*(2*x - 3)*(2*x+1)) /((2*x-1)**3);
            return j_g;
        }

        for (let i = 0;i<phi_range.length;i++){
            data_eq.push(equation_j_g(phi_range[i]));
        }
        
        for (let i=0;i<data_eq.length;i++){
           data_root.push(find_phi(g*data_eq[i],g,root_prec)); 
        }

        let eq_line = {
            x: phi_range,
            y: data_eq,
            type: 'scatter',
            name: 'J curve equation',
        };
        let root_line = {
            x: data_root,
            y: data_eq,
            type: 'scatter',
            name: 'J curve root',
        };
        plot_data.push(eq_line);
        plot_data.push(root_line);
        
        //console.log("J");
        //console.log(plot_data);
        //console.log(plot_data.length);
        return plot_data;
    }
/*
    function produce_plot_colour(){

        $("#j-display").html($("input#jnorm").val());//update display

        let j = parseFloat($("input#jnorm").val());

        let g = 1e3
        let rootprec = 10**(-8);
        let h = 0.1;

        let plot_data = [];

        let data_forwards = results(h,j,g,rootprec);  //(n,h,j,g,root_prec)
        let data_backwards = results(-h,j,g,rootprec);  //(n,h,j,g,root_prec)

        rho_data_plot_b = data_backwards[0];
        phi_data_plot_b = data_backwards[1];

        rho_data_plot_b.reverse();
        phi_data_plot_b.reverse();
        rho_data_plot_b.pop();
        phi_data_plot_b.pop();

        rho_data_plot_f = data_forwards[0];
        phi_data_plot_f = data_forwards[1];
        
        let rho_data = rho_data_plot_b.concat(rho_data_plot_f);
        let phi_data = phi_data_plot_b.concat(phi_data_plot_f);

        //console.log("rho_data");
        rho_data = math.divide(rho_data,h);
        //rho_data = math.round(rho_data,1)// number of decimal places of h
        //console.log(rho_data);

        let length_array = rho_data.length;

        //let theta_array = numeric.linspace(0,2*Math.pi(),rho_data.length);
        let phi_array_quarter = make2Darray(length_array,length_array);

        for(let i = 0; i< length_array; i++){
            for(let j = i; j < length_array; j++){
                //console.log("in loop");
                let r = ((rho_data[i])**2 + (rho_data[j])**2)**0.5;
                //console.log(r);
                //let r_round = math.round(r);
                let rho_val = math.round(r);//number of decimal places of h
                //console.log(rho_val);
                //let rho_index = rho_data.indexOf(rho_val);
                //console.log(rho_index);
                let phi_val = phi_data[rho_val];
                phi_array_quarter[i][j] = phi_val;
                phi_array_quarter[j][i] = phi_val;
            }
        }
        //console.log("quarter array");
        //console.log(phi_array_quarter);
        //let phi_array = make2Darray(2*length_array,2*length_array);
        let data = {
            z: phi_array_quarter,
            type: 'contour'
            };
        plot_data.push(data);
        //console.log(data);
        return plot_data;
    }

    function make2Darray(cols,rows){
        let arr = new Array(cols);
        for (let i = 0; i< arr.length; i++){
            arr[i] = new Array(rows);
        }
        //console.log("arr");
        //console.log(arr);
        return arr;
    };
*/
    function produce_plot_floating(){

        let plot_data = [];
        let g;
        
        //let j_range = [0.1,0.3,1,3,10];
        //let j_range = numeric.linspace(0.1,20,100);
        let j_range = logspace(-4,4,200);
        //console.log("hi");
        //console.log(j_range);

        let rootprec = 10**(-8);
        let h;
        let rho_up_lim;
        let rho_lower_lim = 0;
        let Z = 1;
        let alpha = (Z**0.5)*(1836/(4*Math.PI));
        //console.log(alpha);
        //let min_diff = 0.005;
        let P_array = [];
        let n_s_phi_array = [];

        for (let i = 0; i< j_range.length; i++){
            //console.log(i);
            //console.log(j_range[i]);
            if (j_range[i] < 1){
                g = 1e3;
                h = 0.01;
                rho_up_lim = 15;

            }
            else{
                g = 1e3;
                h = 0.01;
                rho_up_lim = 15;
            }


            let data_forwards = results(h,j_range[i],g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)#
            let data_backwards = results(-h,j_range[i],g,rootprec,rho_up_lim,rho_lower_lim);  //(n,h,j,g,root_prec)
            rho_data_plot_b = data_backwards[0];
            phi_data_plot_b = data_backwards[1];

            rho_data_plot_b.reverse();
            phi_data_plot_b.reverse();
            rho_data_plot_b.pop();
            phi_data_plot_b.pop();

            rho_data_plot_f = data_forwards[0];
            phi_data_plot_f = data_forwards[1];
            //console.log(rho_data_plot_b);
            //console.log(rho_data_plot_f);
            

            let data = [rho_data_plot_b.concat(rho_data_plot_f), phi_data_plot_b.concat(phi_data_plot_f)];
            //console.log(data);
            /*
            let v_line = {
                x: data[0],
                y: data[1],
                type: 'lines',
                name: 'Phi curve J = ' + j_range[i].toString(),
                line: {
                    dash: 'dashdot',
                    width: 4
                  }
            };
            plot_data.push(v_line);
            */

            
            let n_s_phi = math.log(    math.divide(     math.multiply(      math.dotMultiply(data[0],data[0]),    alpha),      j_range[i])      );
            /*
            let v_s_line = {
                x: data[0],
                y: n_s_phi,
                type: 'lines',
                name: 'Phi Surface curve J = ' + j_range[i].toString(),
                line: {
                    dash: 'dot',
                    width: 4
                  }
            };
           
            plot_data.push(v_s_line);
            */

             
            //console.log("found;" + i.toString());
            let diff_NaN = math.abs(math.subtract(data[1],n_s_phi));

            //console.log(diff_NaN);
            let diff_no_NaN  = diff_NaN.filter(function (value) {
                return !Number.isNaN(value);
            });
            //console.log("min");
            //console.log(Math.min.apply(diff));
            //console.log(diff_no_NaN);
            let diff = math.abs(diff_no_NaN);

            /*
            let found = diff.find(function(element){
                return element < min_diff;
                });
            */
                
            let found = Math.min(...diff);//new spread operator
            //console.log(found);
            let rho_index = diff.indexOf(found);
            //console.log(rho_index);

            let P_val = data[0][rho_index];

            P_array.push(P_val);
            n_s_phi_array.push(n_s_phi[rho_index]);
            };

        //console.log("plot data");
        //console.log(P_array);
        //console.log(n_s_phi_array);
    
        let n_s_line = {
            x: P_array,
            y: n_s_phi_array,
            type: 'scatter',
            name: 'Phi Surface curve',
        };
        plot_data.push(n_s_line);
        return plot_data;
    }

    function update_graph(){//update animation
        //console.log("update");
        //console.log();

                //console.log("phi");
        Plotly.animate("graph_phi",
            {data: produce_plot()},//updated data
            {
                fromcurrent: true,
                transition: {duration: 0,},
                frame: {duration: 0, redraw: false,},
                //mode: "afterall"
                mode: "immediate"
            }
        );
        /*
        Plotly.animate("graph_phi_colour",
            {data: produce_plot_colour()},//updated data
            {
                fromcurrent: true,
                transition: {duration: 0,},
                frame: {duration: 0, redraw: false,},
                //mode: "afterall"
                mode: "immediate"
            }
        ); */
        Plotly.animate("graph_floating",
            {data: produce_plot_floating()},//updated data
            {
                fromcurrent: true,
                transition: {duration: 0,},
                frame: {duration: 0, redraw: false,},
                //mode: "afterall"
                mode: "immediate"
            }
        ); 
          
    }

    function initial(){
        
        Plotly.purge("graph_J");
        Plotly.newPlot("graph_J", produce_J_plot(),plt.layoutJ);
        //console.log("j plot");

        Plotly.purge("graph_phi");
        Plotly.newPlot("graph_phi", produce_plot(),plt.layoutKA);
        //console.log("phi plot");

        Plotly.purge("graph_floating");
        Plotly.newPlot("graph_floating", produce_plot_floating(),plt.layoutFloating);

        //dom.gSlider.on("change", update_graph());
        //dom.jSlider.on("change", update_graph());
        //dom.dSlider.on("change", update_graph());

    }
    initial();
});