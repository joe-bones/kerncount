// kerncount launcher: kerncount.exe (Windows) and kerncount.app (Mac, Apple silicon). No admin rights needed.
//
// First start (needs internet, a few minutes, ~1 GB):
//  1. gets uv (a single-file Python manager) into the user folder
//  2. uv fetches Python 3.12 and creates an environment
//     Windows: %LOCALAPPDATA%\kerncount + %LOCALAPPDATA%\kce (short path: TensorFlow's files break
//     Windows' 260-character limit when nested deep). Mac: ~/Library/Application Support/kerncount
//  3. installs the kerncount source embedded in this launcher, plus its dependencies
//
// Every start: opens the kerncount web interface in the browser. Closing the console/Terminal window stops kerncount.
// On the Mac, a double-click from Finder has no window, so the launcher re-opens itself in Terminal.
// A new kerncount.exe (different embedded source) updates the environment automatically.
//
// Flags: --reinstall (rebuild the environment), any other arguments go to `kerncount gui`
// (e.g. a workspace folder). If the folder containing kerncount.exe has a data\ folder, it is the workspace.
package main

import (
	"archive/tar"
	"archive/zip"
	"bufio"
	"bytes"
	"compress/gzip"
	"crypto/sha256"
	_ "embed"
	"encoding/hex"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"time"
)

//go:embed app.zip
var appZip []byte

const uvBase = "https://github.com/astral-sh/uv/releases/latest/download/"

func uvURL() string {
	if runtime.GOOS == "darwin" {
		return uvBase + "uv-aarch64-apple-darwin.tar.gz"
	}
	return uvBase + "uv-x86_64-pc-windows-msvc.zip"
}

func main() {
	code := run()
	if code == -1 { // handed over to a Terminal window
		os.Exit(0)
	}
	if code != 0 {
		fmt.Println("\nPress Enter to close this window.")
		bufio.NewReader(os.Stdin).ReadString('\n')
	}
	os.Exit(code)
}

func run() int {
	exe, _ := os.Executable()
	exe, _ = filepath.EvalSymlinks(exe)
	if runtime.GOOS == "darwin" && !isTerminal() {
		// started from Finder: show progress and keep a window to close, as on Windows
		if err := exec.Command("/usr/bin/open", "-a", "Terminal", exe).Run(); err == nil {
			return -1
		}
	}
	fmt.Println("kerncount")
	fmt.Println("=========")
	base, envDir := dirs()
	exeDir := filepath.Dir(exe)
	// inside kerncount.app/Contents/MacOS: the folder that holds the .app counts
	if app := filepath.Dir(filepath.Dir(exeDir)); strings.HasSuffix(app, ".app") {
		exeDir = filepath.Dir(app)
	}

	var guiArgs []string
	reinstall := false
	for _, a := range os.Args[1:] {
		if a == "--reinstall" {
			reinstall = true
		} else {
			guiArgs = append(guiArgs, a)
		}
	}

	sum := sha256.Sum256(appZip)
	version := hex.EncodeToString(sum[:])[:12]
	py := pythonPath(envDir)
	stamp := filepath.Join(envDir, "kerncount-app.version")
	have, _ := os.ReadFile(stamp)
	if reinstall || strings.TrimSpace(string(have)) != version || !exists(py) {
		if err := install(base, envDir, version, reinstall); err != nil {
			fmt.Println("\nInstall failed:", err)
			fmt.Println("Behind a university/company proxy? Check you can open https://github.com in a browser, then run kerncount.exe again.")
			fmt.Println("To start completely fresh, run the launcher with --reinstall")
			return 1
		}
		os.WriteFile(stamp, []byte(version), 0o644)
	}

	fmt.Println("Starting the kerncount web interface...")
	cmd := exec.Command(py, append([]string{"-m", "kerncount.cli", "gui"}, guiArgs...)...)
	cmd.Dir = exeDir // a data\ folder next to kerncount.exe is picked up as the workspace
	cmd.Stdin, cmd.Stdout, cmd.Stderr = os.Stdin, os.Stdout, os.Stderr
	cmd.Env = append(os.Environ(), "PYTHONIOENCODING=utf-8")
	if err := cmd.Run(); err != nil {
		if ee, ok := err.(*exec.ExitError); ok {
			return ee.ExitCode()
		}
		fmt.Println("Could not start kerncount:", err)
		return 1
	}
	return 0
}

func dirs() (base, env string) {
	home, _ := os.UserHomeDir()
	switch runtime.GOOS {
	case "windows":
		local := os.Getenv("LOCALAPPDATA")
		if local == "" {
			local = filepath.Join(home, "AppData", "Local")
		}
		return filepath.Join(local, "kerncount"), filepath.Join(local, "kce")
	case "darwin":
		base = filepath.Join(home, "Library", "Application Support", "kerncount")
		return base, filepath.Join(base, "env")
	default:
		base = filepath.Join(home, ".local", "share", "kerncount")
		return base, filepath.Join(base, "env")
	}
}

// isTerminal: started from a terminal window? Finder gives /dev/null as stdin, which is also a character
// device, and no TERM variable.
func isTerminal() bool {
	if os.Getenv("TERM") == "" {
		return false
	}
	fi, err := os.Stdin.Stat()
	if err != nil || fi.Mode()&os.ModeCharDevice == 0 {
		return false
	}
	if null, err := os.Stat(os.DevNull); err == nil && os.SameFile(fi, null) {
		return false
	}
	return true
}

func pythonPath(env string) string {
	if runtime.GOOS == "windows" {
		return filepath.Join(env, "Scripts", "python.exe")
	}
	return filepath.Join(env, "bin", "python")
}

func exists(p string) bool { _, err := os.Stat(p); return err == nil }

func install(base, envDir, version string, fresh bool) error {
	fmt.Println("Setting up kerncount (first start or update). This needs internet and takes a few minutes.")
	fmt.Println("Nothing is installed system-wide; everything goes into", base, "and", envDir)
	if err := os.MkdirAll(base, 0o755); err != nil {
		return err
	}
	src := filepath.Join(base, "app-"+version)
	os.RemoveAll(src)
	if err := unzip(appZip, src); err != nil {
		return fmt.Errorf("unpacking kerncount: %w", err)
	}
	// remove sources of older versions
	old, _ := filepath.Glob(filepath.Join(base, "app-*"))
	for _, o := range old {
		if o != src {
			os.RemoveAll(o)
		}
	}
	uv, err := findUV(base)
	if err != nil {
		return err
	}
	if fresh {
		fmt.Println("Removing the old environment...")
		os.RemoveAll(envDir)
	}
	steps := [][]string{
		{"venv", "--python", "3.12", "--seed", "--allow-existing", envDir},
		{"pip", "install", "--python", pythonPath(envDir), src},
	}
	for i, s := range steps {
		fmt.Printf("\n[%d/%d] uv %s\n", i+1, len(steps), strings.Join(s[:2], " "))
		if err := runUV(uv, s); err != nil {
			return err
		}
	}
	fmt.Println("\n[check] loading TensorFlow and StarDist...")
	chk := exec.Command(pythonPath(envDir), "-c", "import tensorflow, stardist, kerncount; print('kerncount OK, TensorFlow', tensorflow.__version__)")
	chk.Stdout, chk.Stderr = os.Stdout, os.Stderr
	if err := chk.Run(); err != nil {
		fmt.Println("\nWARNING: TensorFlow does not load. kerncount still opens, and the 'classical' detector works")
		fmt.Println("(Settings -> Finding nuclei -> Detector).")
		if runtime.GOOS == "windows" {
			fmt.Println("For the recommended StarDist detector, this PC needs the Microsoft Visual C++ Redistributable")
			fmt.Println("(ask IT): https://aka.ms/vs/17/release/vc_redist.x64.exe")
		}
		time.Sleep(4 * time.Second)
	}
	return nil
}

func runUV(uv string, args []string) error {
	cmd := exec.Command(uv, args...)
	cmd.Stdout, cmd.Stderr = os.Stdout, os.Stderr
	env := os.Environ()
	if os.Getenv("UV_NATIVE_TLS") == "" {
		env = append(env, "UV_NATIVE_TLS=1") // trust the Windows certificate store (university proxies)
	}
	cmd.Env = env
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("uv %s failed: %w", args[0], err)
	}
	return nil
}

// findUV: next to kerncount.exe, on PATH, a user install, or our own copy; else download it.
func findUV(base string) (string, error) {
	name := "uv"
	if runtime.GOOS == "windows" {
		name = "uv.exe"
	}
	exe, _ := os.Executable()
	home, _ := os.UserHomeDir()
	cands := []string{filepath.Join(filepath.Dir(exe), name), filepath.Join(base, name), filepath.Join(home, ".local", "bin", name)}
	if p, err := exec.LookPath(name); err == nil {
		cands = append([]string{p}, cands...)
	}
	for _, c := range cands {
		if exists(c) {
			return c, nil
		}
	}
	if runtime.GOOS != "windows" && runtime.GOOS != "darwin" {
		return "", fmt.Errorf("uv not found: install it from https://docs.astral.sh/uv/")
	}
	dst := filepath.Join(base, name)
	fmt.Println("Downloading uv (one file, ~20 MB)...")
	data, err := download(uvURL())
	if err != nil {
		// Go only uses proxy environment variables; PowerShell / curl can do more.
		fmt.Println("Direct download failed (", err, "); trying another way...")
		tmp := filepath.Join(base, "uv-download")
		var dl *exec.Cmd
		if runtime.GOOS == "windows" {
			dl = exec.Command("powershell", "-NoProfile", "-Command",
				fmt.Sprintf("[Net.ServicePointManager]::SecurityProtocol='Tls12'; Invoke-WebRequest -UseBasicParsing -Uri '%s' -OutFile '%s'", uvURL(), tmp))
		} else {
			dl = exec.Command("curl", "-fL", "-o", tmp, uvURL())
		}
		dl.Stdout, dl.Stderr = os.Stdout, os.Stderr
		if perr := dl.Run(); perr != nil {
			return "", fmt.Errorf("could not download uv. Download %s yourself, unpack it, and put %s next to kerncount", uvURL(), name)
		}
		data, err = os.ReadFile(tmp)
		os.Remove(tmp)
		if err != nil {
			return "", err
		}
	}
	if strings.HasSuffix(uvURL(), ".tar.gz") {
		return dst, untarFile(data, name, dst)
	}
	zr, err := zip.NewReader(bytes.NewReader(data), int64(len(data)))
	if err != nil {
		return "", fmt.Errorf("uv download is not a zip: %w", err)
	}
	for _, f := range zr.File {
		if filepath.Base(f.Name) == name {
			if err := writeZipFile(f, dst); err != nil {
				return "", err
			}
			return dst, nil
		}
	}
	return "", fmt.Errorf("uv.exe not found in the download")
}

func untarFile(data []byte, name, dst string) error {
	gz, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return fmt.Errorf("uv download is not a .tar.gz: %w", err)
	}
	tr := tar.NewReader(gz)
	for {
		h, err := tr.Next()
		if err == io.EOF {
			return fmt.Errorf("%s not found in the download", name)
		}
		if err != nil {
			return err
		}
		if filepath.Base(h.Name) == name && h.Typeflag == tar.TypeReg {
			w, err := os.OpenFile(dst, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o755)
			if err != nil {
				return err
			}
			if _, err := io.Copy(w, tr); err != nil {
				w.Close()
				return err
			}
			return w.Close()
		}
	}
}

func download(url string) ([]byte, error) {
	client := &http.Client{Timeout: 5 * time.Minute}
	resp, err := client.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()
	if resp.StatusCode != 200 {
		return nil, fmt.Errorf("HTTP %d", resp.StatusCode)
	}
	return io.ReadAll(resp.Body)
}

func unzip(data []byte, dst string) error {
	zr, err := zip.NewReader(bytes.NewReader(data), int64(len(data)))
	if err != nil {
		return err
	}
	root := filepath.Clean(dst) + string(os.PathSeparator)
	for _, f := range zr.File {
		p := filepath.Join(dst, filepath.FromSlash(f.Name))
		if !strings.HasPrefix(p, root) {
			return fmt.Errorf("bad path in archive: %s", f.Name)
		}
		if f.FileInfo().IsDir() {
			os.MkdirAll(p, 0o755)
			continue
		}
		if err := writeZipFile(f, p); err != nil {
			return err
		}
	}
	return nil
}

func writeZipFile(f *zip.File, dst string) error {
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return err
	}
	r, err := f.Open()
	if err != nil {
		return err
	}
	defer r.Close()
	w, err := os.OpenFile(dst, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o755)
	if err != nil {
		return err
	}
	if _, err := io.Copy(w, r); err != nil {
		w.Close()
		return err
	}
	return w.Close()
}
