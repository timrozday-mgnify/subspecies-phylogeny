process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data' :
        'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc' }"

    input:
    path versions, stageAs: "?/*"

    output:
    path "software_versions.yml",     emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 << 'PYEOF'
import yaml, glob

tools = {}
for f in sorted(glob.glob('**/*.yml', recursive=True)):
    try:
        with open(f) as fh:
            data = yaml.safe_load(fh)
        if isinstance(data, dict):
            for proc, vals in data.items():
                if isinstance(vals, dict):
                    for tool, ver in vals.items():
                        tools[tool] = str(ver)
    except Exception:
        pass

with open('software_versions.yml', 'w') as f:
    yaml.dump({'Software Versions': tools}, f, default_flow_style=False)

items = ''.join(
    f'<dt>{k}</dt><dd><samp>{v}</samp></dd>'
    for k, v in sorted(tools.items())
)
mqc = {
    'id': 'software_versions',
    'section_name': 'Pipeline Software Versions',
    'plot_type': 'html',
    'description': 'Collected at run time from the software output.',
    'data': '<dl class="dl-horizontal">' + items + '</dl>',
}
with open('software_versions_mqc.yml', 'w') as f:
    yaml.dump(mqc, f, default_flow_style=False)
PYEOF
    """

    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml
    """
}
