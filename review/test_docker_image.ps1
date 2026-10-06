param([string]$Image = 'teoloup/hematology-aml-flt3-itd:asv-fixes-20260923')
$ErrorActionPreference = 'Stop'
$repoPath = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$resultPath = Join-Path $repoPath 'review/runs/docker_20260923'
New-Item -ItemType Directory -Force -Path $resultPath | Out-Null
function Invoke-CheckedDocker([string[]]$DockerArgs, [string]$Log) {
    ConvertTo-Json -InputObject $DockerArgs -Compress | Add-Content -Encoding UTF8 (Join-Path $resultPath 'docker_commands.jsonl')
    # Windows PowerShell wraps native stderr (including INFO logs) as errors.
    $ErrorActionPreference = 'Continue'
    & docker @DockerArgs *> $Log
    $dockerExit = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    if ($dockerExit -ne 0) { throw "Docker failed; see $Log" }
}
Invoke-CheckedDocker -DockerArgs @('run','--rm','--entrypoint','python',
    '-e','OMP_NUM_THREADS=1','-e','OPENBLAS_NUM_THREADS=1',
    '-v',"${repoPath}/tests:/tests:ro",$Image,
    '-m','unittest','discover','-s','/tests','-v') -Log (Join-Path $resultPath 'unit_tests.log')
Invoke-CheckedDocker -DockerArgs @('run','--rm','--entrypoint','Nano_ITDseeker',$Image,'--help') -Log (Join-Path $resultPath 'help.log')
foreach ($case in @(
    @{ Name='sim_no_wt'; Folder='review/synthetic' },
    @{ Name='10808_hg38_RG'; Folder='bam_data/test_bam' },
    @{ Name='14417_2runs_hg38_RG'; Folder='bam_data/test_bam' }
)) {
    $sample = $case.Name
    $inputPath = Join-Path $repoPath $case.Folder
    Invoke-CheckedDocker -DockerArgs @('run','--rm',
        '-e','OMP_NUM_THREADS=1','-e','OPENBLAS_NUM_THREADS=1',
        '-v',"${inputPath}:/input:ro",'-v',"${resultPath}:/output",$Image,
        '-b',"/input/$sample.bam",'-o',"/output/$sample",'-s',$sample,
        '-g','hg38','-t','4','--min-allele-frequency','0.02','--html-report') -Log (Join-Path $resultPath "$sample.log")
}
