/*jshint esversion: 6 */
$(window).on('load', function() {//main  
    const dom = {
        gSlider: $("input#gamma"),//gamma slider
        jSlider: $("input#j"),//j slider
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
            layoutColour: {
                autosize: true,
                margin: {
                    l: 40, r: 40, b: 40, t: 40, pad: 10
                },
            }
        };

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

    function results(h,j,g,root_prec){
        //console.log(j);
        let current_step = inital_conditions(j,g,root_prec);
        //console.log("current_step - [inital_rho,inital_u,inital_v]");
        //console.log(current_step);
        let phi_data = [];
        let rho_data = [];

        if (h > 0 && current_step[0] >= 15){
            return [phi_data,rho_data];
        }

        rho_data.push(current_step[0]);
        phi_data.push(current_step[1]);
        
        if(h > 0){
            while(current_step[0] > Math.abs(h) && current_step[0] <= 15){
                current_step = step(current_step[0],current_step[1],current_step[2],h,j);//(rho,u,v,h,J)
                rho_data.push(current_step[0]);
                phi_data.push(current_step[1]);
                
            } 
        }else{
            while(current_step[0] > Math.abs(h)){
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
        let h = 0.001;

        for (let i = 0; i< j_range.length; i++){
            //console.log(j_range[i]);
            let data_forwards = results(h,j_range[i],g,rootprec);  //(n,h,j,g,root_prec)
            let data_backwards = results(-h,j_range[i],g,rootprec);  //(n,h,j,g,root_prec)
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

    function produce_plot_colour(){
        $("#j-display").html($("input#j").val());//update display

        let j = parseFloat($("input#j").val());

        let g = 1e3
        let rootprec = 10**(-8);
        let h = 0.01;

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
        //let theta_array = numeric.linspace(0,2*Math.pi(),rho_data.length);
        phi_array = [];
        
        console.log("shifting array");
        console.log(rho_data.length);
        
        for(let i = 0; i< rho_data.length; i++){
            phi_array.push([]);
            for(let j = 0; j< rho_data.length; j++){
                let r = math.round((i**2 + j**2)**0.5);
                phi_array[i][j] = phi_data[r];
            }   
        }

        let data = [{
            z: phi_array,
            type: 'contour'
            }
        ];
        return data;
    }

    function update_graph(identifier){//update animation
        switch(identifier) {
            case "phi":
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
              // code block
                break;
            case "phi_colour":
                Plotly.animate("graph_phi_colour",
                    {data: produce_plot_colour()},//updated data
                    {
                        fromcurrent: true,
                        transition: {duration: 0,},
                        frame: {duration: 0, redraw: false,},
                        //mode: "afterall"
                        mode: "immediate"
                    }
                ); 
                break;             // code block
          }
    }

    function initial(){
        
        Plotly.purge("graph_J");
        Plotly.newPlot("graph_J", produce_J_plot(),plt.layoutJ);
        console.log("j plot");

        Plotly.purge("graph_phi");
        Plotly.newPlot("graph_phi", produce_plot(),plt.layoutKA);
        console.log("phi plot");

        Plotly.purge("graph_phi_colour");
        Plotly.newPlot("graph_phi_colour", produce_plot_colour(),plt.layoutColour);
        console.log("phi colour plot");

        dom.gSlider.on("change", update_graph("phi"));
        dom.jSlider.on("change", update_graph("phi_colour"));

    }
    initial();
});