/*jshint esversion: 6*/
$(window).on('load', function() {//main
    
    function plot_data(){
        //returns data
    }
    function update_graph(){//update animation
            
        Plotly.animate("graph",
            {data: plot_data()},//updated data
            {
                fromcurrent: true,
                transition: {duration: 0,},
                frame: {duration: 0, redraw: false,},
                //mode: "afterall"
                mode: "immediate"
            }
        );
    }

    function play_loop(){//adds time evolution
        if(isPlay === true) {
            t++;
            Plotly.animate("graph",
                {data: plot_data()},
                {
                    fromcurrent: true,
                    transition: {duration: 0,},
                    frame: {duration: 0, redraw: false,},
                    //mode: "afterall"
                    mode: "immediate"
                });
            requestAnimationFrame(play_loop);//loads next frame
        }
        return 0;
    }

    function initial() {
        Plotly.purge("graph");
        Plotly.newPlot('graph', plot_data(),plt.layout);//create animation

        //dom.tswitch.on("change", update_graph);
        //dom.aSlider.on("input", update_graph);
        //dom.afSlider.on("input", update_graph);

        $('#playButton').on('click', function() {
            document.getElementById("playButton").value = (isPlay) ? "Play" : "Stop";//change play/stop label
            isPlay = !isPlay;
            t = 0;//reset time
            requestAnimationFrame(play_loop);
        });
    }
initial();
});