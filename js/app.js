// Virüs konfigürasyonları
const configurations = {
    sars2: {
        name: "SARS-CoV-2",
        pdb: "data/7K8S.pdb",
        fasta: "data/sars2_sequences.fasta",
        top3Json: "data/sars2_final_top3.json",
        weights: "$$S_{final} = 0.62 \\cdot Stats + 0.38 \\cdot AI$$",
        structuralDesc: "<strong>7K8S (Spike Proteini)</strong> 3D analizi. Renkli bölgeler, tespit edilen farklı hedeflerdir.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    },
    sars1: {
        name: "SARS-CoV-1",
        pdb: "data/5X58.pdb",
        fasta: "data/sars1_sequences.fasta",
        top3Json: "data/sars1_final_top3.json",
        weights: "$$S_{final} = 0.50 \\cdot Stats + 0.50 \\cdot AI$$",
        structuralDesc: "<strong>5X58</strong> 3D analizi. SARS-CoV-1 için evrimsel korunmuş bölgeler.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    },
    flu: {
        name: "Influenza A",
        pdb: "data/1RVX.pdb",
        fasta: "data/flu_sequences.fasta",
        top3Json: "data/flu_final_top3.json",
        weights: "$$S_{final} = 0.55 \\cdot Stats + 0.45 \\cdot AI$$",
        structuralDesc: "<strong>1RVX (Hemagglutinin)</strong> 3D analizi. Antijenik sürüklenmeye dirençli adaylar.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    }
};

const candidateColorsRGB = [
    {r: 239, g: 68, b: 68},  // Kırmızı
    {r: 59,  g: 130, b: 246}, // Mavi
    {r: 16,  g: 185, b: 129}  // Yeşil
];

let viewerInstance = null;
let currentTop3Data = [];
let currentVirusKey = 'sars2';

$(function () {
    viewerInstance = new PDBeMolstarPlugin();
    updateInterface('sars2');
});

async function updateInterface(virusKey) {
    currentVirusKey = virusKey;
    const config = configurations[virusKey];
    document.getElementById('formula-display').innerText = config.weights;
    document.getElementById('structural-desc').innerHTML = config.structuralDesc;
    if (window.MathJax) MathJax.typesetPromise();

    const grid = document.getElementById('top-candidates-grid');
    grid.innerHTML = '<p>Yükleniyor...</p>';

    try {
        // Sekans dosyasını oku
        const fastaResponse = await fetch(config.fasta);
        const fastaText = await fastaResponse.text();
        let fullSequence = '';
        fastaText.split('\n').forEach(line => {
            line = line.trim();
            if (line && !line.startsWith('>')) fullSequence += line;
        });

        // JSON verisini oku
        const response = await fetch(config.top3Json);
        currentTop3Data = await response.json();
        
        grid.innerHTML = ''; 
        currentTop3Data.forEach((candidate, index) => {
            const rgb = candidateColorsRGB[index];
            const hexColor = `rgb(${rgb.r}, ${rgb.g}, ${rgb.b})`;
            const card = document.createElement('div');
            card.className = 'metric-card';
            card.style.borderTopColor = hexColor; 

            // Dizilimi doğrudan JSON'dan al
            const seqStr = candidate.seq_snippet || "N-A";

            // PDB residue id'si doğrudan pos'tan geliyor
            const startPosPdb = candidate.pos - 2;
            const endPosPdb = candidate.pos + 2;

            card.innerHTML = `
                <span class="badge" style="background-color: ${hexColor}; color: white; border: none;">Aday ${index+1}</span>
                <h4>${config.typeMapping[index]}</h4>
                <div style="background-color: #f8fafc; padding: 0.5rem; border-radius: 4px; border: 1px solid #e2e8f0; font-family: monospace; font-size: 1.1rem; text-align: center; margin: 10px 0; color: ${hexColor}; font-weight: bold; letter-spacing: 2px;">
                    ${seqStr}
                </div>
                <p style="font-size: 0.9rem; margin-bottom: 0.5rem;">Pozisyon Aralığı: <strong>${startPosPdb} - ${endPosPdb}</strong> (Merkez: ${candidate.pos})</p>
                <p style="font-size: 0.9rem;">Erişilebilirlik (SASA): ${candidate.sasa} Å²<br>AI Skoru: %${(candidate.ai_score * 100).toFixed(1)}</p>
                <div class="metric-score" style="margin-top: 1rem;">${candidate.final_score.toFixed(3)}<span>S<sub>final</sub></span></div>
            `;
            grid.appendChild(card);
        });
    } catch (err) { grid.innerHTML = `<p>Hata: ${err.message}</p>`; }

    loadProteinStructure(config.pdb);
}

function loadProteinStructure(pdbPath) {
    const options = {
        customData: { url: pdbPath, format: 'pdb' },
        assemblyId: '1', // Biyolojik asambleyi (tekil trimer) göstermeye zorlar
        bgColor: {r: 255, g: 255, b: 255},
        hideControls: false,
        hideSequencePanel: false, // Molstar sekans panelini göster
        sequencePanel: true,
        landscape: true,
        visualStyle: 'cartoon'
    };
    
    const viewerContainer = document.getElementById('myViewer');
    viewerInstance.render(viewerContainer, options);
    
    // Sekans panelini açık tutmaya zorla
    setTimeout(() => {
        if (viewerInstance && viewerInstance.plugin && viewerInstance.plugin.layout) {
            viewerInstance.plugin.layout.setProps({ 
                showSequence: true,
                regionState: { top: 'show', right: 'show' } 
            });
        }
        renderMolstarStyles();
    }, 2000);

    // Yedek event
    if (viewerInstance.events && viewerInstance.events.loadComplete) {
        viewerInstance.events.loadComplete.subscribe(() => {
            setTimeout(renderMolstarStyles, 1000);
        });
    }
}

function renderMolstarStyles() {
    if (!viewerInstance || !viewerInstance.visual) return;
    
    let selectionData = [];
    
    currentTop3Data.forEach((candidate, index) => {
        for (let i = -2; i <= 2; i++) {
            selectionData.push({
                struct_asym_id: 'A',
                residue_number: candidate.pos + i,
                color: candidateColorsRGB[index],
                focus: false
            });
        }
    });

    viewerInstance.visual.select({
        data: selectionData,
        nonSelectedColor: {r: 220, g: 220, b: 220}
    });
}

function startAutomatedAnalysis() {
    const selected = document.getElementById('virus-select').value;
    const overlay = document.getElementById('loading-overlay');
    overlay.style.display = 'flex';
    let step = 0;
    const interval = setInterval(() => {
        step++;
        if (step >= 5) {
            clearInterval(interval);
            overlay.style.display = 'none';
            updateInterface(selected);
        }
    }, 500);
}
