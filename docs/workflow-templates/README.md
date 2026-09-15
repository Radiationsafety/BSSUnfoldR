# Pending GitHub Actions workflows

These three workflow files were prepared locally but **could not be pushed**
because the Personal Access Token used to create the repository does not have
the `workflow` scope (GitHub refuses to push `.github/workflows/*.yaml` via a
PAT without that scope).

## How to enable CI / pkgdown / README rendering

### Option A — Update the PAT (recommended)

1. Open https://github.com/settings/tokens
2. Click the **Generate new token (classic)** button (or edit the existing one).
3. Tick at minimum: `repo` + `workflow`.
4. Generate / save.
5. Locally:

   ```bash
   cd <your-clone>
   git mv docs/workflow-templates/R-CMD-check.yaml   .github/workflows/
   git mv docs/workflow-templates/pkgdown.yaml        .github/workflows/
   git mv docs/workflow-templates/render-readme.yaml  .github/workflows/
   git commit -m "Enable CI / pkgdown / render-README workflows"
   git push origin main
   ```

### Option B — Create the workflow files via the GitHub web UI

1. Open https://github.com/Radiationsafety/BSSUnfoldR/new/main/.github/workflows
2. For each of the three templates below, create a file with the same name
   and copy-paste the contents from `docs/workflow-templates/`.
3. Commit directly to `main` via the UI (you are the repo owner, so this is
   allowed without a PAT).

After either option, also enable **Settings → Pages → Build and deployment →
Source = "GitHub Actions"** so the pkgdown site can publish.

## Workflow summary

| File                       | Purpose                                              | Triggers on              |
| -------------------------- | ---------------------------------------------------- | ------------------------ |
| `R-CMD-check.yaml`         | Run `R CMD check` on Ubuntu/macOS/Windows (R release + devel) | push to main, PRs |
| `pkgdown.yaml`              | Build and deploy the pkgdown documentation site      | push to main, releases, manual |
| `render-readme.yaml`        | Re-render `README.Rmd` → `README.md`                 | changes to README.Rmd    |
