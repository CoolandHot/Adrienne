const myChart = echarts.init(document.getElementById("container"), "dark");
const smallChart = echarts.init(document.getElementById('small_container'), "dark");
var graph_data, separated_json, merged_json;

separated_json = document.querySelector("#detail-pathway").value
merged_json = document.querySelector("#merge-node-pathway").value

var chart_option = {
    title: {
        text: "Pathway connections",
    },
    toolbox: {
        show: true,
        feature: {
            myTool1: {
                show: true,
                title: 'upload json data',
                icon: 'path://M512 0c-282.77 0-512 229.23-512 512s229.23 512 512 512 512-229.23 512-512-229.23-512-512-512zM768 576h-256v192h-128v-192h-256l384-384 384 384z',
                onclick: function () {
                    var input = document.createElement('input');
                    input.type = 'file';
                    input.onchange = function () {
                        console.log("load " + input.files[0].name)
                        make_graph(input)
                    };
                    input.click();
                }
            },
            dataView: { readOnly: true },
            saveAsImage: {}
        }
    },
    tooltip: {
        formatter: function (params) {
            return `${params.name}: ${params.value}<br>${'proportion' in params.data ? params.data.proportion : ''}`
        }
    },
};




// *******************************
// ******** custom functions *****
// *******************************
/**
 * when clicking on the edge of the mergeNode graph, 
 * display a smaller graph for both end nodes
 * @param {*} json_data 
 * @param {*} dom the small popup window for displaying the chart
 */
const update_small_chart = (json_data) => {
    const option = {
        // title: { text: json_data.title, textStyle: { fontSize: 12 } },
        tooltip: {},
        legend: { data: json_data.categories.map(function (a) { return a.name; }), bottom: "5%" },
        series: [
            {
                type: 'graph',
                layout: 'force',
                draggable: true,
                roam: true,
                data: json_data.nodes,
                links: json_data.links,
                categories: json_data.categories
            }
        ]
    };
    smallChart.setOption(option)
}


const update_dataList_input = () => {
    const dataList = document.getElementById("pathway_names")
    graph_data.nodes.map(o => {
        const newOption = document.createElement('option');
        newOption.value = o['name'];
        dataList.appendChild(newOption);
    })
}

var make_graph = async (json_path) => {
    // read the json file
    if (json_path instanceof HTMLInputElement) {
        const file = json_path.files[0];
        graph_data = await new Promise((resolve) => {
            const reader = new FileReader();
            reader.onload = () => {
                resolve(JSON.parse(reader.result));
            };
            reader.readAsText(file);
        });
    } else {
        const response = await fetch(json_path);
        graph_data = await response.json();
    }

    // set new title & legend
    if (graph_data.hasOwnProperty('title')) {
        chart_option['title'] = {
            text: graph_data.title,
        }
    }
    chart_option['legend'] = {
        data: graph_data.categories.map((a) => {
            return a.name;
        }),
        bottom: "3%",
        left: "1%"
    }

    // data series
    chart_option['series'] = [
        {
            type: "graph",
            layout: 'force',
            draggable: true,
            categories: graph_data.categories,
            nodes: graph_data.nodes,
            links: graph_data.links.filter(l => l.value > slider.value),
            emphasis: {
                focus: 'adjacency',
                label: {
                    position: 'right',
                    show: true
                },
                edgeLabel: {
                    show: true,
                    formatter: "{c} common genes"
                }
            },
            roam: true,
            lineStyle: {
                width: 0.5,
                curveness: 0.3,
                opacity: 0.7
            }
        },
    ]

    // update graph
    myChart.setOption(chart_option);

    // update the filter input datalist
    update_dataList_input()
}

const update_edge_threshold = (shared_gene_threshold) => {
    myChart.clear();
    myChart.showLoading();

    chart_option['series'][0]['links'] = graph_data.links.filter(l => l.value > shared_gene_threshold)

    myChart.setOption(chart_option);
    myChart.hideLoading();
}

make_graph(merged_json)


// *******************************
// ********** listeners **********
// *******************************
const slider = document.querySelector('#threshold_slider');
slider.addEventListener('input', () => {
    var threshold = slider.value;
    document.querySelector("label[for='threshold_slider']").innerText = threshold;
    update_edge_threshold(threshold)
});

// json files
document.querySelector("#detail-pathway").addEventListener("change", (e) => {
    separated_json = e.target.value
})
document.querySelector("#merge-node-pathway").addEventListener("change", (e) => {
    merged_json = e.target.value
    make_graph(merged_json)
})

// highlight the filtered input
document.getElementById("filter_pathway").addEventListener('change', (e) => {
    let input = e.target.value.trim()
    let found_pathways = graph_data.nodes.filter(o => o.name.includes(input));
    if (found_pathways.length > 0) {
        const match_indices = found_pathways.map(o => graph_data.nodes.indexOf(o));
        myChart.dispatchAction({
            type: 'highlight',
            seriesIndex: 0,
            dataIndex: match_indices
        });
    }
})

// toggle between the two json files
document.getElementById("nodeMerge_switch").addEventListener("change", (e) => {
    if (e.target.checked) make_graph(merged_json)
    else make_graph(separated_json)
})

// hide the popup panel when clicking outside of it
document.getElementById("panel").addEventListener('click', function (event) {
    if (!document.querySelector('.panel-content').contains(event.target)) {
        document.getElementById("panel").classList.add("hidden");
    }
});

// click on the node/edge to show another detail graph
myChart.on('click', async (params) => {
    // https://echarts.apache.org/handbook/en/concepts/event
    if (params.dataType === 'node') {
        document.getElementById("panel").classList.remove("hidden");

        const response = await fetch(separated_json);
        const json_data = await response.json();
        let new_links = json_data.links.filter(o => o['source'].includes(params.name) | o['target'].includes(params.name))
        let filtered_json = {
            categories: json_data.categories,
            links: new_links,
            nodes: json_data.nodes.filter(o => new_links.some(n => o['name'].includes(n['source']) | o['name'].includes(n['target']))),
        }
        update_small_chart(filtered_json)
        document.querySelector(".panel-content h3").innerText = params.name
    }
    if (params.dataType === 'edge') {
        document.getElementById("panel").classList.remove("hidden");
        const response = await fetch(separated_json);
        const json_data = await response.json();
        const end_nodes = params.name.split(" > ");
        let filtered_json = {
            categories: json_data.categories,
            nodes: json_data.nodes.filter(o => end_nodes.some((n) => o['name'].includes(n))),
            links: json_data.links.filter(o => end_nodes.every((n) => o['source'].includes(n) | o['target'].includes(n)))
        }
        update_small_chart(filtered_json)
        document.querySelector(".panel-content h3").innerText = `${params.name.replace(">", "<===>")}(${params.value} common genes)`
    }
});

window.addEventListener('resize', function () {
    myChart.resize();
});