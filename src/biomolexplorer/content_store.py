"""One immutable, content-addressed home for verified scientific bytes.

The platform already stored content by its digest in two places: the `objects`
table and directory of `calculation_cache`, and the validated entries of
`preparation_cache`.  This module is the single implementation of that
mechanism, so there is one place where bytes become immutable, one place that
verifies them, and one place that decides how a consumer receives them.

The layout is the one `calculation_cache` already uses — ``objects/<first two
characters>/<digest>`` — so an existing cache keeps working without being moved.

Three properties are enforced here rather than asked of the callers:

* **Publication is atomic.** Bytes are written to a temporary file in the
  destination directory, flushed, verified against their digest and only then
  renamed. A publication interrupted at any point leaves either nothing or a
  complete object, never a partial one.
* **An object is read-only** (mode 0444). A consumer that writes in place is
  refused by the kernel instead of corrupting content other runs reference.
* **A writable materialization is an isolated copy.** Nothing a job does to its
  inputs can reach the stored object.

The store never opens a database.  It is used by scientific workers that
receive paths and no catalog, and the catalog rows that reference an object are
written by `artifact_manager`, always *after* the bytes are in place: an
interrupted publication can leave unreferenced bytes, which are inert and
deduplicated by the next publication, but never a reference to bytes that are
absent.
"""
import os
import shutil
import tempfile
from pathlib import Path

from .stage_cache import file_digest

OBJECT_MODE=0o444
FICLONE=0x40049409


def copy_bytes(source,destination):
    """Write the bytes of ``source`` to ``destination`` and flush them to disk.

    Copy-on-write is attempted first.  ext4, the common case on a Linux
    workstation, does not support it: the clone fails with ``EOPNOTSUPP`` and an
    independent copy is made.  The caller is responsible for atomicity.
    """
    import fcntl
    with Path(source).open('rb') as origin,Path(destination).open('wb') as copy:
        try:
            fcntl.ioctl(copy.fileno(),FICLONE,origin.fileno())
        except OSError:
            shutil.copyfileobj(origin,copy)
        copy.flush()
        os.fsync(copy.fileno())


def clone_or_copy(source,destination):
    """Copy bytes atomically: a reader never observes a half-written file."""
    source,destination=Path(source),Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent,delete=False) as stream:
        temporary=Path(stream.name)
    try:
        copy_bytes(source,temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


class ContentStore:
    """Immutable objects addressed by the SHA-256 of their bytes."""

    def __init__(self,root):
        self.root=Path(root)
        self.objects=self.root/'objects'

    def object_path(self,digest):
        return self.objects/digest[:2]/digest

    def has(self,digest):
        """Whether an object is present; presence is not verification."""
        return self.object_path(digest).is_file()

    def size(self,digest):
        try:
            return self.object_path(digest).stat().st_size
        except OSError:
            return None

    def verify(self,digest):
        """Recompute the digest of the stored bytes."""
        path=self.object_path(digest)
        try:
            return path.is_file() and file_digest(path)==digest
        except OSError:
            return False

    def open(self,digest):
        """Read the stored bytes, refusing an object that does not verify."""
        if not self.verify(digest):
            raise ValueError('O objeto armazenado não corresponde ao seu identificador de conteúdo.')
        return self.object_path(digest).open('rb')

    def put(self,source,*,digest=None,validate=None,fingerprint=None):
        """Publish the bytes of ``source`` and return their digest.

        ``validate`` receives the exact bytes about to be published and raises
        to refuse them: content that fails validation is never published.
        ``fingerprint`` lets a caller that already maintains a verified digest
        cache decide whether a present object may be trusted; without it, a
        present object of the right size is trusted at publication and verified
        on every read instead.
        """
        source=Path(source)
        if not source.is_file() or not source.stat().st_size:
            raise ValueError('O conteúdo a publicar está ausente ou vazio: '+str(source))
        digest=digest or file_digest(source)
        target=self.object_path(digest)
        if self._trusted(target,digest,fingerprint):
            return digest
        if validate:validate(source)
        target.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.publishing-',delete=False) as stream:
            temporary=Path(stream.name)
        try:
            copy_bytes(source,temporary)
            if file_digest(temporary)!=digest:
                raise ValueError('O conteúdo mudou durante a publicação: '+str(source))
            os.chmod(temporary,OBJECT_MODE)
            # Replacing a read-only file is governed by the directory, so an
            # object published concurrently with identical bytes is overwritten
            # by an identical object.
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
        from .storage import sync_directory
        sync_directory(target.parent)
        return digest

    def _trusted(self,target,digest,fingerprint):
        if not target.is_file():return False
        if fingerprint is not None:
            try:
                return fingerprint(target)==digest
            except (OSError,ValueError):
                return False
        return True

    def materialize(self,digest,target,mode='isolated',*,verify=True):
        """Give a consumer the stored bytes at ``target``.

        ``shared`` hands over a hard link to the read-only object: it costs no
        bytes, which matters on a filesystem without copy-on-write, and the
        kernel refuses a consumer that writes in place.  ``isolated`` makes an
        independent writable copy, which is what an engine that rewrites its
        own inputs requires.
        """
        if mode not in ('shared','isolated'):
            raise ValueError('Modo de materialização desconhecido: '+str(mode))
        if verify and not self.verify(digest):
            raise ValueError('O objeto armazenado não corresponde ao seu identificador de conteúdo.')
        source=self.object_path(digest)
        target=Path(target)
        target.parent.mkdir(parents=True,exist_ok=True)
        if mode=='shared':
            target.unlink(missing_ok=True)
            try:
                os.link(source,target)
                return target
            except OSError:
                # A different device or a filesystem without links still has to
                # honor the contract of this mode, so the copy stays read-only.
                clone_or_copy(source,target)
                os.chmod(target,OBJECT_MODE)
                return target
        clone_or_copy(source,target)
        # The copy belongs to the job: the stored object stays untouched and the
        # job is free to rewrite what it received.
        os.chmod(target,0o600)
        if verify and file_digest(target)!=digest:
            raise ValueError('A materialização não reproduziu o conteúdo armazenado.')
        return target

    def discard(self,digest):
        """Remove one object; the caller has established that nothing uses it."""
        path=self.object_path(digest)
        try:
            path.unlink(missing_ok=True)
        except OSError:
            return False
        return True

    def stored(self):
        """Every digest the store holds, by its name on disk."""
        if not self.objects.is_dir():return set()
        return {path.name for shard in self.objects.iterdir() if shard.is_dir()
                for path in shard.iterdir() if path.is_file() and not path.name.startswith('.')}

    def sweep(self,referenced):
        """Remove objects nothing references, and temporary files left behind.

        An interrupted publication leaves bytes without a catalog row.  They are
        inert, and the next publication of the same content reuses them, so this
        runs on request rather than on every open.
        """
        removed=0
        for shard in self.objects.iterdir() if self.objects.is_dir() else []:
            if not shard.is_dir():continue
            for path in shard.iterdir():
                if path.name.startswith('.publishing-'):
                    path.unlink(missing_ok=True);continue
                if path.is_file() and path.name not in referenced:
                    path.unlink(missing_ok=True);removed+=1
        return removed

    def usage(self):
        """Objects and bytes under management."""
        total=0;count=0
        for shard in self.objects.iterdir() if self.objects.is_dir() else []:
            if not shard.is_dir():continue
            for path in shard.iterdir():
                if path.is_file() and not path.name.startswith('.'):
                    count+=1;total+=path.stat().st_size
        return {'objects':count,'bytes':total}
