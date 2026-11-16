# One-click kerncount: `kerncount.exe` (Windows) and `kerncount.app` (Mac, Apple silicon)

No admin rights, no Python, no terminal commands. Put the launcher anywhere (ideally in your workspace folder,
next to `data`) and double-click it.

- **First start** (needs internet, a few minutes, ~1 GB): downloads uv, Python 3.12 and kerncount's
  dependencies into your user folder. Nothing system-wide.
  - Windows: `%LOCALAPPDATA%\kerncount` and `%LOCALAPPDATA%\kce`
  - Mac: `~/Library/Application Support/kerncount`
- **Every start**: the kerncount page opens in your browser. A console window (Windows) or Terminal window (Mac)
  shows progress; keep it open while you work and close it (or press **Quit** on the page) to stop.
- **Workspace**: if the launcher sits next to a `data` folder, that's the workspace. Otherwise choose or create one
  on the page.
- **Update**: replace the launcher with a newer build; the environment updates itself on the next start.
- **Repair**: start it with `--reinstall` from a terminal
  (`kerncount.exe --reinstall`, or `kerncount.app/Contents/MacOS/kerncount --reinstall`).

## First-time security prompts (the launchers are not signed with a paid developer certificate)

- **Windows**: "Windows protected your PC" → **More info → Run anyway**.
- **Mac**: unzip `kerncount-mac.zip`, then double-click `kerncount.app`. macOS says it can't verify the developer:
  click **Done**, open **System Settings → Privacy & Security**, scroll down and click **Open Anyway** next to
  "kerncount", then confirm. Only needed once. (Alternative, in Terminal: `xattr -dr com.apple.quarantine kerncount.app`.)
  Allowing it to control Terminal, if asked, is needed for the progress window.

If TensorFlow will not load (Windows: missing Visual C++ runtime), kerncount still opens and the `classical`
detector works.

## Building them

Needs Go ≥ 1.22 on any OS (cross-compiles): `sh launcher/build.sh` →
`dist/kerncount.exe`, `dist/kerncount.app` and `dist/kerncount-mac.zip` (share the zip, not the bare .app).
Or on GitHub: Actions → *Build kerncount launchers* → Run workflow, then download the artifacts.
The launchers embed the current `src/` at build time, so rebuild after changing the code.
