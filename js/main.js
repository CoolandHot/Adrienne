var file_paras = [
    { "file": "", "name": "", "reverse": true, "up_regulate": 'down' },
    { "file": "", "name": "", "reverse": true, "up_regulate": 'down' },
    { "file": "", "name": "", "reverse": true, "up_regulate": 'down' },
    { "file": "", "name": "", "reverse": true, "up_regulate": 'down' },
    { "file": "", "name": "", "reverse": true, "up_regulate": 'down' },
];
var p_adjust = 0.05;


function readCSVFile(file, reverse = true) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = event => {
            const csvData = event.target.result;
            const rows = csvData.split(/\r?\n/);
            const headers = rows.shift().replaceAll('"', '').split(',');
            headers[0] = 'gene'; // the first header is empty
            const records = rows.map(row => {
                const values = row
                    .replaceAll('"', '')
                    .split(',')
                const record = {};
                headers.forEach((header, index) => {
                    if (reverse & (header == "log2FoldChange" | header == "stat")) {
                        record[header] = -parseFloat(values[index]);
                    } else {
                        record[header] = values[index];
                    }
                });
                return record;
            });

            const filteredRecords = records.filter(record => record['padj'] < p_adjust);
            resolve(filteredRecords);
        };
        reader.onerror = error => reject(error);
        reader.readAsText(file);
    });
}



// ===========================================================
// ************************** echarts ***********************
// ===========================================================

var d3Container = d3.select("#container");
var vennChart = venn.VennDiagram()
    .width(500)
    .height(500);
var tooltip = d3.select("body")
    .append("div")
    .attr("class", "venntooltip");

const map_csv = () => {
    return new Promise((resolve, reject) => {
        var myPromises = file_paras
            .filter(v => v['file'] != "")
            .map(async v => {
                let records = await readCSVFile(v['file'], v['reverse']);
                // only reserve the up/down-regulation
                records = records.filter(record => v['up_regulate'] == 'up' ? record['log2FoldChange'] > 0 : record['log2FoldChange'] < 0);
                // get the remaining gene list
                let geneList = records.map(record => record['gene']);
                return { "name": `${v['name']}(${v['up_regulate']})`, "value": geneList, "records": records }
            });
        Promise.all(myPromises).then((results) => {
            resolve(results);
        });
    });

}

// *********************************
// ******* make intersection *******
// *********************************
function intersectList(...arrays) {
    return arrays.reduce((a, b) => a.filter(c => b.includes(c)));
}

const compute_intersect = (dataset) => {

    // singleton
    let data_one = dataset.map(v => ({ sets: [v['name']], size: v['value'].length, records: v['records'] })), data_two = [], data_three = [], data_four = [], data_five = [];

    // intersection of two sets
    let data_two01 = intersectList(dataset[0]['value'], dataset[1]['value']), data_two02, data_two12, data_two03, data_two13, data_two23;
    data_two = [
        {
            sets: [dataset[0]['name'], dataset[1]['name']],
            size: data_two01.length,
            records: data_two01
        }
    ];

    // intersection of three/four sets
    if (dataset.length > 2) {
        data_two02 = intersectList(dataset[0]['value'], dataset[2]['value']);
        data_two12 = intersectList(dataset[1]['value'], dataset[2]['value']);
        data_two = data_two.concat([
            { sets: [dataset[0]['name'], dataset[2]['name']], size: data_two02.length, records: data_two02 },
            { sets: [dataset[1]['name'], dataset[2]['name']], size: data_two12.length, records: data_two12 }
        ])
        let data_three012 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[2]['value']);
        data_three = data_three.concat([
            {
                sets: dataset.slice(0, 3).map(v => v['name']),
                size: data_three012.length,
                records: data_three012
            }
        ])
    }
    if (dataset.length > 3) {
        data_two03 = intersectList(dataset[0]['value'], dataset[3]['value']);
        data_two13 = intersectList(dataset[1]['value'], dataset[3]['value']);
        data_two23 = intersectList(dataset[2]['value'], dataset[3]['value']);
        data_two = data_two.concat([
            { sets: [dataset[0]['name'], dataset[3]['name']], size: data_two03.length, records: data_two03 },
            { sets: [dataset[1]['name'], dataset[3]['name']], size: data_two13.length, records: data_two13 },
            { sets: [dataset[2]['name'], dataset[3]['name']], size: data_two23.length, records: data_two23 }
        ])
        let data_three013 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[3]['value']);
        let data_three023 = intersectList(dataset[0]['value'], dataset[2]['value'], dataset[3]['value']);
        let data_three123 = intersectList(dataset[1]['value'], dataset[2]['value'], dataset[3]['value']);
        data_three = data_three.concat([
            {
                sets: [0, 1, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three013.length,
                records: data_three013
            },
            {
                sets: [0, 2, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three023.length,
                records: data_three023
            },
            {
                sets: [1, 2, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three123.length,
                records: data_three123
            },
        ])

        let data_four0123 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[2]['value'], dataset[3]['value']);
        data_four = data_four.concat([
            {
                sets: dataset.map(v => v['name']),
                size: data_four0123.length,
                records: data_four0123
            }
        ])
    }
    if (dataset.length > 4) {
        let data_two04 = intersectList(dataset[0]['value'], dataset[4]['value']);
        let data_two14 = intersectList(dataset[1]['value'], dataset[4]['value']);
        let data_two24 = intersectList(dataset[2]['value'], dataset[4]['value']);
        let data_two34 = intersectList(dataset[3]['value'], dataset[4]['value']);
        data_two = data_two.concat([
            { sets: [dataset[0]['name'], dataset[4]['name']], size: data_two04.length, records: data_two04 },
            { sets: [dataset[1]['name'], dataset[4]['name']], size: data_two14.length, records: data_two14 },
            { sets: [dataset[2]['name'], dataset[4]['name']], size: data_two24.length, records: data_two24 },
            { sets: [dataset[3]['name'], dataset[4]['name']], size: data_two34.length, records: data_two34 }
        ])

        let data_three014 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[4]['value']);
        let data_three024 = intersectList(dataset[0]['value'], dataset[2]['value'], dataset[4]['value']);
        let data_three034 = intersectList(dataset[0]['value'], dataset[3]['value'], dataset[4]['value']);
        let data_three124 = intersectList(dataset[1]['value'], dataset[2]['value'], dataset[4]['value']);
        let data_three134 = intersectList(dataset[1]['value'], dataset[3]['value'], dataset[4]['value']);
        let data_three234 = intersectList(dataset[2]['value'], dataset[3]['value'], dataset[4]['value']);
        data_three = data_three.concat([
            {
                sets: [0, 1, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three014.length,
                records: data_three014
            },
            {
                sets: [0, 2, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three024.length,
                records: data_three024
            },
            {
                sets: [0, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three034.length,
                records: data_three034
            },
            {
                sets: [1, 2, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three124.length,
                records: data_three124
            },
            {
                sets: [1, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three134.length,
                records: data_three134
            },
            {
                sets: [2, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three234.length,
                records: data_three234
            },
        ])

        let data_four0124 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[2]['value'], dataset[4]['value']);
        let data_four0134 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[3]['value'], dataset[4]['value']);
        let data_four0234 = intersectList(dataset[0]['value'], dataset[2]['value'], dataset[3]['value'], dataset[4]['value']);
        let data_four1234 = intersectList(dataset[1]['value'], dataset[2]['value'], dataset[3]['value'], dataset[4]['value']);

        data_four = data_four.concat([
            {
                sets: [0, 1, 2, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_four0124.length,
                records: data_four0124
            },
            {
                sets: [0, 1, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_four0134.length,
                records: data_four0134
            },
            {
                sets: [0, 2, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_four0234.length,
                records: data_four0234
            },
            {
                sets: [1, 2, 3, 4].map(ind => dataset[ind]).map(v => v['name']),
                size: data_four1234.length,
                records: data_four1234
            },
        ])

        let data_five01234 = intersectList(dataset[0]['value'], dataset[1]['value'], dataset[2]['value'], dataset[3]['value'], dataset[4]['value']);
        data_five = data_five.concat([
            {
                sets: dataset.map(v => v['name']),
                size: data_five01234.length,
                records: data_five01234
            }
        ])

    }
    return data_one.concat(data_two).concat(data_three).concat(data_four).concat(data_five)
}

const plot_update = (dataset) => {
    document.getElementById("chart-title").innerText = dataset.map(v => v['name']).join(" - ");

    var data = compute_intersect(dataset);
    d3Container.datum(data).call(vennChart);

    // outline of each circle
    d3Container.selectAll("path")
        .style("stroke-opacity", 0)
        .style("stroke", "#fff")
        .style("stroke-width", 2);

    // add listeners to all the groups to display tooltip on mouseover
    d3Container.selectAll("g")
        .on("mouseover", function (d, i) {
            // sort all the areas relative to the current item
            venn.sortAreas(d3Container, d);
            // Display a tooltip with the current size
            tooltip.transition().duration(100).style("opacity", .9);
            tooltip.text(d.size + " genes");
            // highlight the current path
            var selection = d3.select(this).transition("tooltip").duration(200);
            selection.select("path")
                .style("stroke-width", 3)
                .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                .style("stroke-opacity", 1);
        })

        .on('click', (d, i) => {
            // sort all the areas relative to the current item
            venn.sortAreas(d3Container, d);
            // pop up detailed tables
            document.getElementById("panel").classList.remove("hidden");

            var records = d.records;
            var tableNames = d.sets, table;
            if (tableNames.length == 1) {
                table = '<h3>' + tableNames + '</h3>';
                table += records.length + ' genes<br>';
                table += '<table><thead><tr><th>gene</th><th>log2FoldChange</th><th>p.adjust</th></tr></thead><tbody>';
                for (let j = 0; j < records.length; j++) {
                    table += '<tr>' + `<td>${records[j]['gene']}</td>` + `<td>${records[j]['log2FoldChange']}</td>` + `<td>${records[j]['padj']}</td>` + `</tr>`;
                }
                table += '</tbody></table>'
            } else {
                table = '<h3>' + tableNames.join(" - ") + '</h3>';
                table += records.length + ' genes<br>';
                table += records.join(", ")
            }
            document.querySelector('.panel-content').innerHTML = table;

        })

        .on("mousemove", function () {
            tooltip.style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
        })

        .on("mouseout", function (d, i) {
            tooltip.transition().duration(200).style("opacity", 0);
            var selection = d3.select(this).transition("tooltip").duration(200);
            selection.select("path")
                .style("stroke-width", 0)
                .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                .style("stroke-opacity", 0);
        });

}



// ===========================================================
// ***************** file picker listeners ******************
// ===========================================================
const filePickers = ['filePicker1', 'filePicker2', 'filePicker3', 'filePicker4', 'filePicker5'];
for (const picker of filePickers) {
    document.getElementById(picker).addEventListener('click', () => {
        const input = document.createElement('input');
        input.type = 'file';
        input.accept = '.csv';
        input.onchange = event => {
            const file = event.target.files[0];
            let num = parseInt(picker.replace("filePicker", ""));
            let input_dom = document.getElementById(`file${num}`);
            if (input_dom.value == "") {
                input_dom.value = file.name.replace(".csv", "");
                input_dom.dispatchEvent(new Event("change"));
            }
            file_paras[num - 1]['file'] = file;
        };
        input.click();
    });
}
// ===========================================================
// *********************** input listeners *******************
// ===========================================================

document.querySelectorAll('input[type="text"]')
    .forEach(i => i.addEventListener('change', event => {
        if (event.target.id == "p_adjust") {
            p_adjust = parseFloat(document.getElementById("p_adjust").value)
        }
        if (event.target.id.includes('file')) { // file name
            let num = parseInt(event.target.id.replace("file", ""));
            file_paras[num - 1]["name"] = event.target.value
        }
    }))
document.querySelectorAll('input[type="checkbox"]')
    .forEach(i => i.addEventListener('change', event => {
        file_paras[parseInt(event.target.id.replace("reverse_fc", "")) - 1]["reverse"] = event.target.checked;
    }))
document.querySelectorAll('input[type="radio"]')
    .forEach(i => i.addEventListener('change', event => {
        let num = parseInt(event.target.name.replace("regulation", ""));
        file_paras[num - 1]["up_regulate"] = document.querySelector(`input[name="${event.target.name}"]:checked`).value;
    }))
document.querySelectorAll('.clear-file')
    .forEach(i => i.addEventListener('click', event => {
        let parent = event.target.parentNode;
        let num = parseInt(parent.getAttribute("for").replace("file", ""));
        parent.querySelector("input").value = "";
        file_paras[num - 1]["file"] = "";
        file_paras[num - 1]["name"] = "";

    }))



document.getElementById("panel").addEventListener('click', function (event) {
    if (!document.querySelector('.panel-content').contains(event.target)) {
        document.getElementById("panel").classList.add("hidden");
    }
});

document.getElementById("createChart").addEventListener('click', async event => {
    let data = await map_csv();
    plot_update(data)
});